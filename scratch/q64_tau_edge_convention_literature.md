# Q64 — the τ angular-cell-edge convention: what the literature actually prescribes

Status: **COMPLETE**, plus a **Follow-up** section appended 2026-08-11 that resolves
an internal tension and **CORRECTS two claims in §1/§9**. All three sources read;
every load-bearing equation spot-verified against the rendered scan. Sources are all
local (`scratch/literature/`); no online search was used.

⚠ **Read `## Follow-up` §F0 before acting on §1 or §9.** Two claims there are
refuted by measurement: the "(C) diverges because of a missing weight
renormalisation" hypothesis (the real cause is node/partition incompatibility), and
the "use β as the oracle" recommendation (β is blind on ORPHEUS's own rule, though
it is sound on the sphere).

---

## 0. THE ONE-PARAGRAPH ANSWER

**(C), and it is not close.** Three independent sources — Bailey–Morel–Chang (2010)
Eqs. (12)/(52), and Lathrop (2000) p. 249 — define the angular cell edges by
**accumulating the quadrature weights** in the **radial direction cosine**, seeded
at the level's inward endpoint, and define τ as the **fractional position of the
ordinate's radial cosine inside that cell** (BMC Eqs. 42/74; Lathrop Eq. 23 with
τ = (1+δ)/2). **(A) is nobody's convention** — no source forms an arithmetic mean of
two neighbouring ordinate cosines. **(B) is nobody's convention either** — no source
ever bisects an ω-cell; the only azimuthal half-angle quantity in the literature
(Hébert's `η_{p,q±1/2}`) is defined by a **discrete conservation recursion**, not by
evaluating a trigonometric function at a bisected angle. **Sphere and cylinder use
the SAME convention** (BMC say so in as many words), so ORPHEUS's (C)-sphere /
(A)-cylinder asymmetry is ORPHEUS's own. **No source clamps τ**, and BMC's own
worked S₂ example yields τ₁ ≈ 0.4226 — *below* ORPHEUS's ½ floor. Finally, Lathrop
(2000) p. 245/249 states the exact trade-off that the MMS table is measuring:
**τ = ½ has the BEST local truncation order (O(Δμ²)) and the WRONG diffusion
limit**, while the Morel–Montry τ has worse truncation (O(δΔμ + Δμ²)) and the
correct one. The MMS is measuring the property τ was never chosen to optimise.

---

## 1. VERDICT

**None of (A), (B). The literature prescribes (C) — cumulative-quadrature-weight
cell edges taken in the RADIAL DIRECTION COSINE — for the cylinder, and it is the
SAME convention as the sphere, differing only by a per-level rescaling.**

Bailey, Morel & Chang (2010) NSE 165:2, 149–169 is the paper ORPHEUS cites, and it
treats *both* geometries with one construction:

- **Sphere** (their §IV): cell edges by the recursion `μ_{m+1/2} = μ_{m−1/2} + w_m`,
  seeded `μ_{1/2} = −1`, with `Σ_m w_m = 2` — **Eq. (12), PDF p. 5 / printed p. 152**.
- **Cylinder / R-Z** (their §V): the *same* recursion **per ξ-level**,
  `μ_{m+1/2,n} = μ_{m−1/2,n} + w̄_{m,n}`, seeded `μ_{1/2,n} = −√(1−ξ_n²)`, ending at
  `μ_{M+1/2,n} = +√(1−ξ_n²)`, **with the level's weights renormalised to sum to
  `2√(1−ξ_n²)`** — **Eq. (52), PDF p. 10 / printed p. 157**.
- **τ, both geometries**: the fractional position of the ordinate cosine inside its
  own cell, measured **in the radial cosine μ** — Eq. (42) (sphere, PDF p. 8) and
  Eq. (74) (cylinder, PDF p. 13). Eq. (74) is textually declared "analogous to those
  defined by Eq. (42)".

So the literature's answer to the crux question is: **the asymmetry is ORPHEUS's,
not the literature's.** BMC use one convention for both geometries. What differs
between sphere and cylinder in BMC is only the *interval* the recursion runs over
(`[−1, +1]` vs `[−sinθ_n, +sinθ_n]`) and the corresponding **weight
renormalisation on each level**.

⭐ **The most probable cause of the measured (C) divergence.** The brief describes
ORPHEUS's (C) as "boundaries accumulated from −sinθ" with "cell width proportional
to the ordinate's quadrature weight". BMC's Eq. (52) requires that the per-level
weights be **renormalised so that they sum to `2√(1−ξ_n²) = 2 sinθ_n`**. If the
cylinder's accumulation reuses weights normalised to the sphere's `Σw = 2` (or to
1, or to the full-sphere `4π` convention of their Eq. 51) while starting at
`−sinθ_n`, the recursion **overshoots `+sinθ_n`**, the last cell's width has the
wrong sign or magnitude, and τ leaves `[0,1]` — which is exactly the signature of
the reported `nan`. **(C) as BMC define it is not the same as (C) as measured**
unless the renormalisation is present. This should be re-measured before (C) is
declared refuted.

---

## 2. The equations, in ORPHEUS notation

### 2.0 Notation map

BMC use a *geometry-dependent* symbol set. In their **R-Z** section (their Fig. 2
caption text, PDF p. 9): `ξ = cos θ`, `μ = sin θ cos ω`, `η = sin θ sin ω`.
(The OCR renders the third as `sin θ cos ω`, an obvious OCR duplication — see the
scan-verification note in §6.)

| BMC (R-Z, §V) | meaning | ORPHEUS |
|---|---|---|
| `μ` , `μ_{m,n}` | **radial** direction cosine | `mu_x` = **η** |
| `η` | **azimuthal** direction cosine | `mu_y` = **ξ** |
| `ξ` , `ξ_n` | **axial** direction cosine (= `cos θ`, the level) | `mu_z` |
| `ω` | azimuthal angle within the level | **ω** |
| `n` | ξ-level index (polar level) | polar-level index |
| `m` | ordinate index within the level (the ω-march) | azimuthal index |
| `α_{m±1/2,n}` | angular-redistribution coefficient | `α` |
| `w_{m,n}` , `w̄_{m,n}` | ordinate weight / **level-renormalised** weight | `w` |
| `τ_{m,n}` | weighted-diamond factor | `τ` |

⚠ **Notation collision to hold in mind**: BMC's `μ` in R-Z is ORPHEUS's `η`
(radial), and BMC's `η` in R-Z is ORPHEUS's `ξ` (azimuthal) — the two symbols are
*swapped* relative to the ORPHEUS convention stated in the brief. In BMC's
*spherical* section `μ` is the ordinary 1-D streaming cosine. (This is exactly the
L-001 hazard: BMC's `μ` means the radial cosine in curvilinear contexts.)

### 2.1 Sphere — the cell-edge recursion `[VERIFIED-ON-SCAN p.5]`

`Bailey-Morel-Chang(2010) Eq. (12), PDF p. 5, printed p. 152`

```
η_{m+1/2} = η_{m−1/2} + w_m ,     η_{1/2} = −1 ,     η_{M+1/2} = +1 ,
```
with `Σ_{m=1}^{M} w_m = 2` (their Eq. 11). BMC note that `η_{1/2} = −1` is set
explicitly to start the recursion and `η_{M+1/2} = +1` then follows automatically.

**This is candidate (C), and it is exactly ORPHEUS's spherical arm** (`mu_edge[n+1]
= mu_edge[n] + w[n]`, seeded at −1). **Sub-question 2: CONFIRMED, BMC Eq. (12).**

### 2.2 Cylinder — the cell-edge recursion `[VERIFIED-ON-SCAN p.10]`

`Bailey-Morel-Chang(2010) Eq. (52), PDF p. 10, printed p. 157`

```
η_{m+1/2,n} = η_{m−1/2,n} + w̄_{m,n} ,
η_{1/2,n}   = −√(1 − μ_z,n²)  ( = −sin θ_n ) ,
η_{M+1/2,n} = +√(1 − μ_z,n²)  ( = +sin θ_n ) ,
```
**with the level weights renormalised so that `Σ_m w̄_{m,n} = 2√(1 − μ_z,n²)`.**
BMC's own words: the cell-edge cosines "are obtained as in the spherical-geometry
case except that there is a separate recursion and starting cosine for each
ξ-level" (PDF p. 10).

### 2.3 The weighted-diamond closure — identical in both geometries `[VERIFIED-ON-SCAN pp.5,10]`

Sphere, `Eq. (15), PDF p. 5`; cylinder, `Eq. (53), PDF p. 10`:

```
ψ_m = τ_m ψ_{m+1/2} + (1 − τ_m) ψ_{m−1/2}          (sphere)
ψ_{m,n} = τ_{m,n} ψ_{m+1/2,n} + (1 − τ_{m,n}) ψ_{m−1/2,n}   (cylinder)
```
"where each weighting factor `τ` can take on any value between zero and one. Note
that `τ = 1` gives the step scheme, and `τ = ½` gives the diamond scheme."
(PDF p. 5; repeated verbatim in the cylinder case, PDF p. 10.)

### 2.4 τ — the Morel–Montry weighting factors `[VERIFIED-ON-SCAN pp.8,13]`

Sphere, `Eq. (42), PDF p. 8, printed p. 155`:
```
τ_m = (η_m − η_{m−1/2}) / (η_{m+1/2} − η_{m−1/2})
```
Cylinder, `Eq. (74), PDF p. 13, printed p. 160`:
```
τ_{m,n} = (η_{m,n} − η_{m−1/2,n}) / (η_{m+1/2,n} − η_{m−1/2,n})
```

⭐ **Sub-question 1 answered: the fractional position is taken in the RADIAL
COSINE η (BMC's μ), not in ω and not in the cumulative weight itself.** The
*edges* come from the cumulative weight (§2.1/§2.2); the *fraction* is then a ratio
of cosine differences. These are two separate design decisions and ORPHEUS must
match both.

BMC also state the equivalent form the code cites as its definition —
`Eq. (43), PDF p. 8` (sphere):
```
η_m = τ_m η_{m+1/2} + (1 − τ_m) η_{m−1/2}
```
i.e. **the same weighted-diamond relation applied to the COSINES rather than the
fluxes**, and BMC remark that Eq. (43) is what makes the flux closure Eq. (15)
*exact* when the angular flux is linear in the cosine (their Eq. 1 form).

### 2.5 The condition τ exists to satisfy (why it is not free)

`Eq. (41), PDF p. 8` (sphere) and `Eq. (73)/(75), PDF p. 13` (cylinder):
```
β_sph = ½ Σ_m η_m [ α_{m+1/2} η_{m+1/2} − α_{m−1/2} η_{m−1/2} ] = 0
β_cyl = Σ_n Σ_m η_{m,n} [ α_{m+1/2,n} η_{m+1/2,n} − α_{m−1/2,n} η_{m−1/2,n} ] = 0
```
This β is the Morel–Montry factor; it multiplies a **contamination term in the
second-order current** `J⁽²⁾ = −(1/3σ_t)∇φ⁽¹⁾ + (β/(rσ_t²)) ∂φ⁽⁰⁾/∂r` (Eq. 75).
Setting β = 0 *determines* τ. Consequently **τ is not a tunable accuracy knob — it
is the unique member of the weighted-diamond family that removes the O(ε)
diffusion-limit contamination for an arbitrary standard quadrature set**
(BMC, PDF p. 13: the Morel–Montry scheme "is the only method in the family of
general weighted diamond methods … that forces the β factor to be zero for any
standard quadrature set").

---

## 3. Sphere vs cylinder — SAME convention (BMC)

BMC's cylinder section is written explicitly as the sphere's analogue at every
step: Eq. (52) "as in the spherical-geometry case except …"; Eq. (53) is Eq. (15)
with a level index; Eq. (74) is "analogous to those defined by Eq. (42)"; and their
closing paragraph (PDF p. 13) says the first-order asymptotic results "are
fundamentally the same for the spherical and cylindrical geometry cases", with
β_sph and β_cyl written as one pair in Eq. (75).

⟹ **The (C)-for-sphere / (A)-for-cylinder asymmetry in ORPHEUS has no counterpart
in BMC.** (Hébert and Lathrop still to be checked — §6.)

---

---

## 3b. Hébert (2009) Chapter 3 — a DIFFERENT scheme entirely: τ ≡ ½

⭐ **Hébert never introduces τ, a weighted diamond, or Morel–Montry. His angular
closure is the plain DIAMOND — τ = ½ — in BOTH geometries, and he calls it that.**

| | Hébert's angular closure | Eq. | PDF p. / printed p. | status |
|---|---|---|---|---|
| cylinder | `φ_{p,q,i} = ½(φ_{p,q−1/2,i} + φ_{p,q+1/2,i})` | **(3.406)** | p. 73 / 139 | `[VERIFIED-ON-SCAN p.73]` |
| sphere | `φ_{n,i} = ½(φ_{n−1/2,i} + φ_{n+1/2,i})` | **(3.431)** | p. 77 / 143 | `[OCR-ONLY]` |

Both are printed under the italicised label *"the diamond differencing scheme"*.

⟹ **The two Eqs. the ORPHEUS code cites — (3.437) and (3.439) — are the τ = ½
extrapolations, not a τ definition.** Their angular half is
`φ_{n+1/2,i} = 2 φ_{n,i} − φ_{n−1/2,i}`, which is precisely BMC's Eq. (16)
`ψ_{m+1/2} = (1/τ)ψ_m − ((1−τ)/τ)ψ_{m−1/2}` **evaluated at τ = ½**. (3.437) and
(3.439) differ only in the *spatial* half — inward (`φ_{n,i−1/2} = 2φ_{n,i} −
φ_{n,i+1/2}`, μ<0 sweep) versus outward (`φ_{n,i+1/2} = 2φ_{n,i} − φ_{n,i−1/2}`,
μ>0 sweep). **Neither carries any information about angular cell edges.**

⚠ **Consequence for the measurement**: the row "control: τ ≡ ½ everywhere" in the
brief is not a control — **it is Hébert's actual published scheme**. That it wins on
the MMS is expected and is not evidence against the Morel–Montry τ (see §4).

### Hébert's cylinder does NOT define radial-cosine cell edges at all

Hébert defines only the **azimuthal** half-angle cosines `η_{p,q±1/2}`
(= ORPHEUS `ξ` = `mu_y`), and only because they *are* α. `[VERIFIED-ON-SCAN p.72]`

`Hebert(2009) Eq. (3.396)–(3.399), PDF p. 72, printed p. 138`, attributed in the
text to **Alcouffe and O'Dell** (his ref. [36]):

```
(3.396)   ξ_{p, 2F(p)+1/2} = 0            # ORPHEUS ξ; the zero-weight end point
(3.397)   W_{p,q} η_{p,q} C − W_p ( ξ_{p,q+1/2} C − ξ_{p,q−1/2} C ) = 0
(3.398)   ξ_{p,q+1/2} = ξ_{p,q−1/2} + (W_{p,q} / W_p) · η_{p,q}
(3.399)   α_{p,q±1/2} ≡ W_p · ξ_{p,q±1/2}   ⟹   α_{p,q+1/2} = α_{p,q−1/2} + W_{p,q} η_{p,q}
```
(Hébert's `μ` = ORPHEUS `η` radial; Hébert's `η` = ORPHEUS `ξ` azimuthal;
Hébert's `ξ_p` = ORPHEUS `mu_z` level. The rewrite above is in ORPHEUS symbols.)

His stated derivation of (3.398): the values are chosen **"in such a way as to
preserve the constant flux in a case where an infinite, homogeneous and
non-absorbing medium contains a population of particles, in the absence of
sources"** — i.e. **(3.398) is a DISCRETE CONSERVATION condition, not a geometric
evaluation.** Eq. (3.396) is the only *geometric* input: at the level's end point
the direction lies in the ξ–μ plane so the azimuthal cosine is exactly zero.

Angular range and indexing, `[VERIFIED-ON-SCAN p.72]`: `N/2` axial levels in
`0 ≤ ξ ≤ 1`; each level carries `2F(p)` base points in `0 ≤ ω_{p,q} ≤ π`
(Eq. 3.391). The march's starting direction is the zero-weight point at
**ω = π** (Eq. 3.407, `[VERIFIED-ON-SCAN p.73]`), where the redistribution term
vanishes.

### Hébert's sphere likewise defines α, not edges

`Hebert(2009) Eq. (3.423)–(3.424), PDF p. 76, printed p. 142` `[OCR-ONLY]`:
same constant-flux argument, seeded from **"The first value α is equal to
1 − (−1)² = 0"** (i.e. `α_{1/2} = 1 − μ_{1/2}²` at `μ_{1/2} = −1`), then propagated:
```
α_{n+1/2} = α_{n−1/2} − 2 W_n μ_n
```
The limits `μ_{n±1/2}` appear only as limits of integration in Eq. (3.418); Hébert
**never gives a recursion for their values**, because with τ ≡ ½ they are never
needed. Base points and weights are stated to be **N-point Gauss–Legendre**.

---

## 4. The clamp — NO SOURCE PRESCRIBES, JUSTIFIES, OR MENTIONS ANY LIMITER ON τ

Explicitly: **across Bailey–Morel–Chang (2010) and Hébert (2009) Ch. 3 there is no
clamp, no blend, no positivity fixup, and no lower bound on τ of any kind.** What
they say instead:

1. **The admissible range is `[0, 1]`, not `[½, 1]`.** BMC state it twice, once per
   geometry, in identical words: the weighting factors "can take on any value
   between zero and one", with `τ = 1` the step scheme and `τ = ½` the diamond
   scheme — sphere `[VERIFIED-ON-SCAN p.5]` (text under Eq. 15), cylinder
   `[VERIFIED-ON-SCAN p.10]` (text under Eq. 53).

2. ⛔ **BMC's own worked example produces τ < ½, so a `[½,1]` clamp would corrupt
   it.** For S₂ (Gauss, `μ_{3/2} = 0` forced by Eq. 46), BMC's Eq. (47)
   `[VERIFIED-ON-SCAN p.8]`:
   ```
   τ_1 = η_1 + 1 = 1 − 1/√3 ≈ 0.42265        τ_2 = η_2 = 1/√3 ≈ 0.57735
   ```
   `τ_1 ≈ 0.4226 < ½`. **A `max(0.5, ·)` clamp silently replaces the unique β = 0
   member of the family with something else, i.e. it re-introduces exactly the
   O(ε) diffusion-limit contamination the Morel–Montry τ exists to remove.**
   This is a strong, direct, scan-verified refutation of the lower clamp.

3. **What the literature says about τ = ½ specifically** (BMC p. 8, printed 155,
   `[VERIFIED-ON-SCAN p.8]`, paraphrased): setting every τ to ½ is a *simple
   average* relationship; doing so would force the quadrature set to be the
   midpoint rule (for S₂: `μ_1 = −½, μ_2 = ½`), and BMC say the midpoint rule "is
   not a good choice of quadrature sets in general". Their conclusion is that if
   one wants a **general** quadrature set *and* the correct diffusion equation for
   the first-order flux, one **must** use Eq. (42)'s weighting factors.
   ⟹ **τ = ½ is only self-consistent when the ordinates ARE the midpoints of their
   own cumulative-weight cells.** For Gauss–Legendre they are not.

4. **Accuracy/robustness trade-off, in the literature's own terms**: BMC's Table I
   (sphere, Gauss–Legendre) and **Table II** (R-Z, level-symmetric)
   `[VERIFIED-ON-SCAN p.13]` tabulate β. For R-Z the Morel–Montry weighted diamond
   gives β ≈ round-off (`−4.12E−16` at S₂ … `−1.55E−15` at S₁₆) while plain diamond
   gives `1.42E+00` at S₂, `1.15E−01` at S₄, `2.09E−02` at S₈, `1.12E−02` at S₁₂,
   `8.82E−03` at S₁₆ — decaying but never zero. The *observable* consequence BMC
   attribute to β ≠ 0 is the **flux dip** (unphysical positive slope at the axis /
   origin), or a corresponding unphysical negative slope. They report no
   positivity, stability, or robustness problem with the Morel–Montry τ.

5. **Where the literature DOES discuss a negative-flux fixup, it is SPATIAL, not
   angular.** Hébert's negative-flux condition, Eqs. (3.387)–(3.389)
   (`Δx_i Σ_i < 2|μ_n|`, PDF p. 71 `[OCR-ONLY]`), is about the *spatial* diamond in
   slab geometry. Nothing analogous is applied to the angular variable.

⟹ **The `max(0.5, min(1.0, τ))` clamp in ORPHEUS's cylinder is unattributable to
either source, and item 2 shows it is actively inconsistent with BMC's Eq. (42)/(74).**
If it is load-bearing for stability, that is a *numerical* finding of ORPHEUS's own
and must be documented as such — not cited to BMC or Hébert.

---

## 5. α and τ — do they share one cell partition?

**Short answer: NO, not in the way the brief's hypothesis assumes — and the
hypothesis's premise about α is only half right.**

### 5a. α is `W_p × (azimuthal cosine)`, but that cosine is NOT a geometric evaluation

ORPHEUS's assertion, per the brief: "α_{p,q±1/2} = W_p · (tangential cosine)
evaluated at REAL half-angle boundaries in ω".

- **First half — CORRECT and verbatim Hébert**: `α_{p,q±1/2} ≡ W_p η_{p,q±1/2}`
  is literally the defining sentence under Eq. (3.399) `[VERIFIED-ON-SCAN p.72]`,
  and Hébert's `η` is the azimuthal (tangential) cosine = ORPHEUS `ξ`.
- ⛔ **Second half — REFUTED.** The values `η_{p,q±1/2}` are **not** obtained by
  evaluating `sin θ_p · sin ω_{p,q±1/2}` at a bisected angle. They are generated by
  the recursion (3.398), whose stated justification is preservation of a constant
  flux in an infinite non-absorbing medium. The symbols `ω_{p,q±1/2}` appear in
  Eq. (3.393) only as *labels for where the half-angle FLUX lives*; their numeric
  values are never computed, and no equation in §3.9.3 defines them. The single
  geometric anchor is Eq. (3.396), `ξ_{p,2F(p)+1/2} = 0`, which seeds the recursion.

⟹ **There is no "real ω half-angle boundary" in the source construction at all.**
Candidate **(B)** therefore has no support in either reference: neither BMC nor
Hébert ever bisects the ω-cell.

### 5b. In BMC, α and τ come from TWO SEPARATE recursions over the same nominal cell

This is the sharpest structural fact for the ORPHEUS question. In BMC's cylinder
(and identically in the sphere), there are two accumulations, both seeded and
terminated at the level's end points, and they are **different**:

| object | recursion | what it accumulates | used by | Eq. |
|---|---|---|---|---|
| `α_{m+1/2,n}` | `α_{m+1/2,n} = α_{m−1/2,n} − η_{m,n} w_{m,n}` | the **first moment** `Σ w·η` | the redistribution term | (50) cyl, (11) sph |
| `η_{m+1/2,n}` | `η_{m+1/2,n} = η_{m−1/2,n} + w̄_{m,n}` | the **zeroth moment** `Σ w` | **τ only**, via (74)/(42) | (52) cyl, (12) sph |

So `τ` is *not* derived from `α`'s partition and `α` is *not* derived from `τ`'s.
They are two different exact-integration conditions imposed on the same angular
cell. Answering the brief's conditional directly: **the premise "if α is defined at
ω half-angle boundaries then the same partition defines τ" does not arise**, because
α is not defined at ω half-angle boundaries (§5a) and τ's edges are a separate
weight accumulation (BMC Eq. 52).

### 5c. A consistency identity worth knowing `[MY DERIVATION — later CONFIRMED by Lathrop, see §8a]`

Hébert's sphere seeds `α_{1/2} = 1 − μ_{1/2}²` at `μ_{1/2} = −1`. If one *also*
insisted `α_{n+1/2} = 1 − μ_{n+1/2}²` at every edge, Hébert's (3.424) would give
`μ_{n+1/2}² = μ_{n−1/2}² + 2 W_n μ_n`, whereas BMC's (12) gives
`μ_{n+1/2} = μ_{n−1/2} + W_n`. Equating the two forces
`μ_n = μ_{n−1/2} + W_n/2` — **the ordinate is the midpoint of its weight cell, i.e.
exactly τ = ½.**

⟹ The two constructions coincide **iff τ = ½**, and Gauss–Legendre ordinates are
not weight-cell midpoints, so they genuinely differ. This is the same statement
BMC make in prose on p. 8 (item 3 of §4). It also means **the "α uses one partition,
τ uses another" tension the brief senses is real and is inherent to the method** —
it is not an ORPHEUS bug. What ORPHEUS must not do is *invent a third* partition.

---

## 6. Confidence, verification status, and gaps

### Scan-verification ledger

| claim | eq. | page | status |
|---|---|---|---|
| sphere α recursion + `Σw = 2` | BMC (11) | PDF p. 5 | `[VERIFIED-ON-SCAN p.5]` |
| **sphere cell-edge recursion = cumulative weight from −1** | BMC (12) | PDF p. 5 | `[VERIFIED-ON-SCAN p.5]` |
| sphere weighted-diamond closure; τ∈[0,1]; τ=1 step, τ=½ diamond | BMC (15) | PDF p. 5 | `[VERIFIED-ON-SCAN p.5]` |
| sphere edge-flux extrapolation | BMC (16) | PDF p. 5 | `[VERIFIED-ON-SCAN p.5]` |
| β = 0 condition (sphere) | BMC (41) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| **τ definition (sphere)** | BMC (42) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| **cosine weighted-diamond identity** | BMC (43) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| S₂ worked τ values (`τ₁ = μ₁+1 < ½`) | BMC (47) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| "all τ = ½ forces the midpoint rule" prose | — | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| **cylinder cell-edge recursion + per-level renormalisation to `2√(1−ξ²)`** | BMC (52) | PDF p. 10 | `[VERIFIED-ON-SCAN p.10]` |
| cylinder weighted-diamond closure, τ∈[0,1] | BMC (53) | PDF p. 10 | `[VERIFIED-ON-SCAN p.10]` |
| **τ definition (cylinder)** | BMC (74) | PDF p. 13 | `[VERIFIED-ON-SCAN p.13]` |
| β_sph / β_cyl | BMC (75) | PDF p. 13 | `[VERIFIED-ON-SCAN p.13]` |
| Table II β values (R-Z, level-symmetric) | — | PDF p. 13 | `[VERIFIED-ON-SCAN p.13]` |
| cylinder α: `α = W_p η_half`, conservation recursion, Alcouffe–O'Dell | Héb (3.396)–(3.399) | PDF p. 72 | `[VERIFIED-ON-SCAN p.72]` |
| level/ω indexing, `0 ≤ ω ≤ π`, `2F(p)` points | Héb (3.391) | PDF p. 72 | `[VERIFIED-ON-SCAN p.72]` |
| **cylinder closure is the plain diamond (τ ≡ ½)** | Héb (3.406) | PDF p. 73 | `[VERIFIED-ON-SCAN p.73]` |
| cylinder starting direction at ω = π | Héb (3.407) | PDF p. 73 | `[VERIFIED-ON-SCAN p.73]` |
| sphere α recursion, `α_{1/2} = 1−(−1)²` | Héb (3.423)–(3.424) | PDF p. 76 | `[OCR-ONLY]` |
| **sphere closure is the plain diamond (τ ≡ ½)** | Héb (3.431) | PDF p. 77 | `[VERIFIED-ON-SCAN p.77]` |
| **sphere extrapolations cited by ORPHEUS** | Héb (3.437) | PDF p. 77 | `[VERIFIED-ON-SCAN p.77]` |
| sphere μ>0-sweep average flux | Héb (3.438) | PDF p. 77 | `[VERIFIED-ON-SCAN p.77]` |
| sphere μ<0 starting-direction chain | Héb (3.432)–(3.435) | PDF p. 77 | `[VERIFIED-ON-SCAN p.77]` |
| **δ ↔ weighted diamond; `\|δ\| ≤ 1` iff ordinate in cell** | Lathrop (23),(24) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| **α = 1−μ² only when δ = 0; else + β, `Δβ = −δΔμ²/2`** | Lathrop (22),(25),(26) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| **"only with δ = 0 is the truncation order O(Δμ²)"** | Lathrop after (30) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| Reed–Lathrop cosine selection + diffusion condition | Lathrop (27)–(29) | PDF p. 8 | `[VERIFIED-ON-SCAN p.8]` |
| **cell edges `μ_{m+1/2}=μ_{m−1/2}+Δμ_m`, `ΣΔμ_m = 2`** | Lathrop §IV text | PDF p. 12 | `[VERIFIED-ON-SCAN p.12]` |
| `c ≡ Σ μ_m² Δμ_m`; α closed form | Lathrop (54),(58) | PDF p. 12 | `[VERIFIED-ON-SCAN p.12]` |
| **"δ = 0 … would not have c = ⅔ … not the correct diffusion limit"** | Lathrop after (59) | PDF p. 12→13 | `[VERIFIED-ON-SCAN p.12]` |
| Morel–Montry vs Reed–Lathrop lineage | Lathrop after (59) | PDF p. 12 | `[VERIFIED-ON-SCAN p.12]` |

### Print/OCR slips found

- **BMC Eq. (52), first line, is a TYPO IN THE PRINTED JOURNAL** (not an OCR
  artefact): it reads `μ_{m+1/2,n} = μ_{m+1/2,n} + w̄_{m,n}` on the rendered page.
  The right-hand subscript must be `m−1/2`, as in the sphere's Eq. (12) which is
  printed correctly (`μ_{m+1/2} = μ_{m−1/2} + w_m`, verified p. 5), and as required
  for the recursion to be non-trivial. **Use `m−1/2`.**
- BMC's Fig. 2 caption text (PDF p. 9) is OCR'd as `η = sin θ cos ω`, duplicating
  the μ line. The correct reading is `η = sin θ sin ω` (forced by
  `μ² + η² + ξ² = 1`). Low-stakes, but do not transcribe the OCR form.
- Hébert Eqs. (3.420)–(3.422) carry an OCR index slip (`μ_{m±1/2}` / `W_n` mixed);
  the intended index is `n` throughout. Cosmetic.
- **Lathrop Eq. (30): the sidecar reads `−⅓Δμ[…]`; the rendered page reads
  `−¼Δμ[…]`.** OCR slip. Not load-bearing here (it is the third Taylor coefficient),
  but it is the same `⅓`-vs-fraction OCR failure mode already recorded for
  Adams–Larsen 2002. **Use `−¼`.**

### Confidence

- **HIGH** on: BMC's cylinder τ = Eq. (74) with Eq. (52) edges; BMC using one
  convention for both geometries; ORPHEUS's sphere matching BMC Eq. (12) exactly
  (**two** independent sources: BMC + Lathrop); Hébert prescribing τ ≡ ½ in both
  geometries and defining no radial-cosine edges; **no clamp anywhere in any of the
  three sources**; BMC's S₂ τ₁ < ½; the τ = ½ truncation-vs-diffusion-limit
  trade-off (§7c, two verified statements one page apart).
- **HIGH** on the negative result for (A): no source forms an arithmetic mean of two
  neighbouring ordinate cosines anywhere.
- **HIGH** on the negative result for (B): no source bisects an ω-cell anywhere;
  the only azimuthal half-angle object (Hébert's `η_{p,q±1/2}`) is conservation-
  defined, and Lathrop Eqs. (25)/(26) prove the analogous sphere object departs from
  its geometric value whenever τ ≠ ½.
- **MEDIUM-HIGH** on the cylinder transfer of Lathrop's §7c reasoning: Lathrop
  argues in the sphere only, but BMC's §V shows the R-Z algebra is term-by-term
  analogous (their Eqs. 73/75 vs 40/41), so the trade-off carries. Marked as a
  transfer, not a citation.
- **MEDIUM** on the *cause* of ORPHEUS's (C) divergence being the missing per-level
  weight renormalisation — this is inference from BMC Eq. (52) plus the reported
  `nan`, not a measurement. It is cheap to test and should be tested before (C) is
  written off.

### Gaps — what I could NOT determine

1. **BMC's R-Z is 2-D (r,z); ORPHEUS's is 1-D cylindrical.** Structurally the ω-march
   per level is identical (BMC's ξ-levels are ORPHEUS's `mu_z` levels; BMC's μ is
   ORPHEUS's η; BMC even write the redistribution with the same `α/(r w)` shape).
   But BMC never write the 1-D cylinder. Nothing in their derivation depends on the
   z-term (it is annihilated separately, Eq. 70/72), so the transfer is sound — but
   it is a transfer, and I am flagging it rather than hiding it.
2. **Neither source states the ORDERING/sign convention of the ω-march relative to
   the μ-edge recursion** in enough detail to fix ORPHEUS's index direction by
   citation alone. BMC seed at `μ_{1/2,n} = −√(1−ξ_n²)` (inward-most radial cosine)
   and Hébert seeds at `ω = π` (also inward-most, μ = −sinθ). These agree, but the
   sign of ORPHEUS's α accumulation must be checked against its own redistribution
   sign (BMC Eq. 49 carries `+[…]/(r w)`, Hébert Eq. 3.400 carries `−(1/W)[…]` —
   opposite conventions, hence BMC's `α = α − wμ` vs Hébert's `α = α + wμ`).
3. **Hébert's ref. [36] (Alcouffe & O'Dell)** is the primary source for the
   cylindrical α/η construction. It is NOT in `scratch/literature/`. If the exact
   original wording on the half-angle *cell boundaries* matters, that is the
   document to acquire — it is the only place a statement about ω-cell edges could
   plausibly live. (Hébert presents their result; he does not re-derive it.)
4. **Lathrop (2000) is SPHERE-ONLY** — it cannot adjudicate the cylinder on its own
   authority. Its role here is corroborative (it independently reproduces BMC's
   sphere construction) and explanatory (§7c). Everything cylinder-specific rests
   on BMC §V and Hébert §3.9.3.
5. **Morel & Montry (1984), TTSP 13(5):615** — the *primary* source for the
   weighted-diamond τ, cited by both BMC and Lathrop. **NOT in
   `scratch/literature/`.** BMC and Lathrop both restate its result and neither
   contradicts the other, so nothing here depends on reading it — but if a
   statement about τ's *admissible range* or a limiter is ever going to exist, that
   is the remaining place it could be. Recommend acquiring it to close sub-question
   4 at the primary source. Confirmed present in the reference lists of both local
   papers (Lathrop ref. 9; BMC's Morel–Montry attribution throughout).

### Suggested follow-up reads

- `Hebert(2009)Chapter3.pdf` **Fig. 3.25** (referenced at Eq. 3.396) — the
  zero-weight-point diagram; it is the one figure that shows the level's ω-cell
  layout, and OCR cannot carry it. **This is the highest-value remaining local
  read** if any doubt persists about where the cylinder's angular cells sit.
- `Bailey-Morel-Chang(2010)` **PDF p. 9, Fig. 2** (R-Z coordinate system) and
  **PDF p. 10, Fig. 3** (starting directions, weighted directions, ξ-levels for
  triangular S₆) — the two figures that show the level/march layout.
- Alcouffe & O'Dell (Hébert ref. [36]) and Morel & Montry (1984) — NOT LOCAL; user
  decision to acquire.

---

## 7. Lathrop (2000) — the third source, and the adjudication the MMS cannot give

`Lathrop, "A Comparison of Angular Difference Schemes for One-Dimensional Spherical
Geometry S_N Equations", Nucl. Sci. Eng. 134, 239–264 (2000).`

⚠ **Scope: SPHERICAL ONLY.** The word "cylindrical" does not occur in the paper.
It cannot adjudicate the cylinder directly — but it is the most explicit source in
the corpus on *why* the design space looks the way it does, and every structural
statement it makes is the one BMC then carry over to R-Z verbatim.

### 7a. Lathrop's δ IS τ, affinely: `τ = (1 + δ)/2`

`Eq. (23) and (24), PDF p. 8, printed p. 245` `[VERIFIED-ON-SCAN p.8]`

```
(23)   η_m = ((1+δ)/2) η_{m+1/2} + ((1−δ)/2) η_{m−1/2} = η̄ + δ·Δη/2
(24)   ψ_m = ((1+δ)/2) ψ_{m+1/2} + ((1−δ)/2) ψ_{m−1/2}
```
(Lathrop's `μ` is the 1-D spherical streaming cosine = ORPHEUS `mu_x`; I write `η`
only to keep one symbol per concept across this document. `Δμ_m` is his notation
for the **quadrature weight** — his quadrature set is the pair `{μ_m, Δμ_m}`.)

Eq. (23) is BMC's Eq. (43); Eq. (24) is BMC's Eq. (15). With `τ = (1+δ)/2`:
`δ = 0 ⟺ τ = ½` (diamond) and `δ = 1 ⟺ τ = 1` (step).

⭐ **The bound, in the source's own words**: *"where |δ| ≤ 1 as long as `μ_m` is in
`Δμ`."* `[VERIFIED-ON-SCAN p.8]` ⟹ **τ ∈ [0,1] is a CONSEQUENCE of the ordinate
lying inside its own angular cell — not an imposed limiter.** It is the same range
BMC state. Nothing in any source narrows it to `[½, 1]`.

### 7b. Lathrop independently confirms candidate (C) for the sphere

`PDF p. 12, printed p. 249` `[VERIFIED-ON-SCAN p.12]`, paraphrased: he assumes the
quadrature set `{μ_m, Δμ_m}` is symmetric about μ = 0, that
**`μ_{m+1/2} = μ_{m−1/2} + Δμ_m`** for m = 1…N, that `μ_m` relates to those cell-edge
values by Eq. (23), and that **`Σ_m Δμ_m = 2`, "the length of the μ interval"**.

That is BMC's Eq. (12) written by a different author. **Two independent sources,
same construction.** Combined with BMC Eq. (52), (C) is the corpus's convention.

He also gives the α recursion in closed form, `Eq. (58)` `[VERIFIED-ON-SCAN p.12]`:
```
α_{m+1/2} = −2 Σ_{j=1..m} η_j Δη_j
```
— the integrated form of BMC Eq. (11) / Hébert Eq. (3.424). All three agree.

### 7c. ⭐⭐ THE ADJUDICATION: τ = ½ has the BEST truncation order and the WRONG diffusion limit

This is the single most useful thing in the corpus for the present decision, and
Lathrop states both halves explicitly, one page apart.

**Half 1 — truncation error favours τ = ½.** `PDF p. 8, printed 245`
`[VERIFIED-ON-SCAN p.8]`. After Eq. (30) he concludes, of the weighted-diamond
family, that the third Taylor coefficient is `O(δΔμ + Δμ²)`, the next
`O(δΔμ² + Δμ³)`, "and so on, **so that only with `μ_m = μ̄` (δ = 0) is the truncation
order `O(Δμ²)`**". He adds (Eq. 25/26 discussion) that for Gauss–Legendre sets
**"except near the ends of the μ range, the values of δ are small and decrease with
the order of the quadrature"** — so the ordinary diamond's `O(δΔμ)` error "is as
small as `Δμ²` for most directions".

**Half 2 — the diffusion limit forbids τ = ½ with a general quadrature.**
`PDF p. 12→13, printed 249→250` `[VERIFIED-ON-SCAN p.12]`: *"Note that the standard
angular diamond scheme, δ = 0, implies a cell-centered `μ_m = μ̄`, which if used in
Eq. (18) would not have `c = ⅔` and hence would not give the correct diffusion
limit"* — where `c ≡ Σ_m μ_m² Δμ_m` (Eq. 54) and `c = ⅔` is the **diffusion
condition** (Eq. 29).

⟹ **The two criteria point in opposite directions, and the literature knows it.**
Setting τ = ½ *defines* the ordinates to be weight-cell midpoints; midpoint
ordinates on equal intervals do not satisfy `Σ μ² Δμ = ⅔`, so the discrete system
does not reduce to P₁. Conversely, keeping a good quadrature (Gauss–Legendre,
level-symmetric) forces the ordinates off the cell midpoints, and then only the
Morel–Montry τ of Eq. (42)/(74) restores the correct diffusion limit — at the cost
of `O(δΔμ)` in local truncation.

### 7d. The two historical lineages (do not conflate them)

Lathrop `PDF p. 12, printed 249` `[VERIFIED-ON-SCAN p.12]`, paraphrased:

- **Reed & Lathrop** (his ref. 7): use a weighted diamond *and choose the direction
  cosines* to minimise truncation error (his Eqs. 27/28), which **fixes the
  quadrature** and makes it satisfy the diffusion condition Eq. (29) by
  construction. Drawback he names: that quadrature "does not correctly integrate
  Legendre polynomials higher than P₂", so anisotropic-source moments are wrong.
- **Morel & Montry** (his ref. 9, `Transp. Theory Stat. Phys. 13, 5, 615 (1984)`):
  keep **Gauss–Legendre** (so high-order anisotropic sources integrate correctly),
  and then **derive** `δ_m` from Eq. (23) for that given quadrature. This is the
  lineage BMC formalise, and it is the one ORPHEUS is implementing. It "eliminated
  the central flux dip".

⟹ For ORPHEUS the relevant lineage is Morel–Montry: **the quadrature is given, and
τ is a derived quantity — never a free parameter and never a tuned one.**

### 7e. What Lathrop adds about the `μ = −1` starting direction

`PDF p. 6, printed 243` `[OCR-ONLY]`, paraphrased: when `m = 1` and `μ_{m−1/2} = −1`
the angular-derivative term vanishes, giving a plane-geometry-like equation; he
stresses throughout the desirability of an **initialising computation of the
`μ = −1` angular flux** and of a difference approximation that makes the central
flux the same for all directions. He reports a **36 % centre-flux error** for N = 2
with `μ̄ = −½` when this is mishandled. Same structural role as Hébert's `ω = π`
starting direction (Eq. 3.407) and BMC's `μ_{1/2,n} = −√(1−ξ_n²)` (Eq. 55).

---

## 8. Direct answers to the five sub-questions

| # | question | answer | evidence |
|---|---|---|---|
| **main** | which of (A)/(B)/(C)/(D) for the cylinder? | **(C)** — cumulative-weight edges in the radial cosine, per level, weights renormalised to `2√(1−ξ_n²)` | BMC Eq. (52) `[VERIFIED-ON-SCAN p.10]` |
| **1** | fractional position in η, ω, or weight? | **in η, the radial cosine.** The *edges* come from cumulative weight; the *fraction* is a ratio of η-differences. Two separate decisions. | BMC Eqs. (74)/(42) `[VERIFIED-ON-SCAN pp.13,8]`; Lathrop Eq. (23) `[VERIFIED-ON-SCAN p.8]` |
| **2** | is the sphere really (C)? | **YES — confirmed.** ORPHEUS's `mu_edge[n+1] = mu_edge[n] + w[n]` from −1 is BMC Eq. (12) exactly. Confirmed a second time by Lathrop p. 249. | BMC Eq. (12) `[VERIFIED-ON-SCAN p.5]`; Lathrop `[VERIFIED-ON-SCAN p.12]` |
| **3** | same convention for both, or genuinely different? | **SAME.** BMC derive the cylinder as the sphere's analogue at every step and say the first-order results "are fundamentally the same". **The ORPHEUS asymmetry is ORPHEUS's own.** | BMC Eqs. (52),(53),(74),(75) + prose `[VERIFIED-ON-SCAN pp.10,13]` |
| **4** | any clamping/limiting of τ? | **NO SOURCE MENTIONS ANY LIMITER.** The only bound is τ ∈ [0,1], and it is a *consequence* of the ordinate lying in its cell, not an imposition. **BMC's own S₂ gives τ₁ ≈ 0.4226 < ½**, so a `max(0.5,·)` floor contradicts the source. | BMC Eqs. (15),(53),(47) `[VERIFIED-ON-SCAN pp.5,10,8]`; Lathrop Eq. (23) `[VERIFIED-ON-SCAN p.8]` |
| **5** | does α's partition define τ? | **The premise is refuted twice over.** α is *not* evaluated at a geometric ω half-angle — it is a conservation recursion (Hébert 3.398, Lathrop 19/58). And α equals the geometric edge value `1 − η_{edge}²` **iff δ = 0, i.e. iff τ = ½** (Lathrop Eqs. 22/25/26). So α and τ legitimately use different accumulations over one partition; that is the method, not a bug. | Hébert (3.396)–(3.399) `[VERIFIED-ON-SCAN p.72]`; Lathrop (22),(25),(26) `[VERIFIED-ON-SCAN p.8]` |

### 8a. Lathrop's Eqs. (22)/(25)/(26) — the explicit α-vs-τ statement

This is the literature's own answer to sub-question 5, and it is a *stronger* result
than the derivation I sketched in §5c (which it independently reproduces).
`[VERIFIED-ON-SCAN p.8]`

```
(22)  when η_m = η̄  (i.e. τ = ½):    α_{m+1/2} = 1 − η_{m+1/2}²
(25)  in general:                     α_{m+1/2} = 1 − η_{m+1/2}² + β_{m+1/2}
(26)  with                            β_{m+1/2} − β_{m−1/2} = −δ · Δη²/2 ,  β_{1/2} = 0
```

⟹ **α coincides with the geometric value at the cell edge ONLY in the diamond case.**
Whenever τ ≠ ½ the conservation-derived α departs from `1 − η_edge²` by an amount
that accumulates as `−δΔη²/2` per cell. The cell *partition* is shared; the
*evaluation on it* is not, and cannot be, because α is fixed by discrete
conservation while the edges are fixed by the weights.

**Implication for ORPHEUS**: do NOT "fix" a perceived inconsistency by forcing α and
the τ edges onto a common geometric evaluation. That would silently re-impose δ = 0
(τ = ½) and destroy the conservation property α carries.

---

## 9. What this means for the measurement (read this before acting)

The brief's table is measuring **local truncation error on a manufactured
anisotropic solution**, which is precisely the metric Lathrop shows τ = ½ optimises
(§7c, half 1). It is **blind** to the property τ exists to deliver — the O(ε)
diffusion-limit consistency of the *first-order* flux (BMC β = 0), whose observable
is the **flux dip / anomalous slope at the axis in an optically thick diffusive
problem**. Neither the MMS error norm nor its convergence order can see β.

⟹ **The literature-consistent choice is (C) with the derived, UNCLAMPED τ.** Two
concrete follow-ups fall out, in order of cost:

1. ⭐ **Re-measure (C) with the per-level weight renormalisation of BMC Eq. (52).**
   The reported `nan` is the signature of τ leaving [0,1], which is exactly what
   happens if the accumulation from `−sinθ_n` uses weights that do not sum to
   `2 sinθ_n`. `[HYPOTHESIS — inference from Eq. (52) + the reported nan, not measured]`
   Cheap, and it decides whether (C) was ever really tested.
2. **Adopt the literature's own discriminating diagnostic, not the MMS**: BMC's β
   (their Eq. 75) is a pure function of the quadrature and the α/edge recursions —
   no solve required. Computing `β_cyl` for each candidate convention is a direct,
   cheap, non-MMS oracle: the correct convention gives β at round-off (BMC Table II:
   `−4.12E−16` … `−1.55E−15`), the wrong ones give `O(10⁻²–10⁰)`. **This is the
   adjudicator the MMS cannot be.**
3. The confirmatory physical test is BMC §VI / Lathrop §V: a diffusive sphere or
   cylinder (BMC use `R = 1, σ_t = 0.1, σ_a = 0.05, Q = 1`, reflecting at the axis,
   vacuum elsewhere, scaled by ε), looking at the **sign of the flux slope at the
   origin**, not at an L2 norm.

---

## 10. Notation-conflict register (for whoever writes the theory page)

| symbol | BMC spherical §IV | BMC R-Z §V | Hébert §3.9.3 cyl | Lathrop (sphere) | ORPHEUS |
|---|---|---|---|---|---|
| radial cosine | `μ` | `μ` | `μ` | `μ` | `η` (`mu_x`) |
| azimuthal cosine | — | `η` | `η` | — | `ξ` (`mu_y`) |
| axial cosine / level | — | `ξ`, `ξ_n` | `ξ_p` | — | `mu_z` |
| ordinate weight | `w_m` | `w_{m,n}`, `w̄_{m,n}` | `W_{p,q}` | **`Δμ_m`** | `w` |
| level weight | — | (implicit) | `W_p` | — | — |
| closure factor | `τ_m` | `τ_{m,n}` | *(none; ≡ ½)* | **`δ`, with `τ = (1+δ)/2`** | `τ` |
| redistribution coeff. | `α_{m±1/2}` | `α_{m±1/2,n}` | `α_{p,q±1/2}` | `α_{m±1/2}` | `α` |
| level index | — | `n` | `p` | — | polar |
| in-level index | `m` | `m` | `q` | `m` | azimuthal |

⚠ **The three traps.**
1. **`μ`/`η` are SWAPPED** between the sources and ORPHEUS in curvilinear contexts:
   BMC/Hébert `μ` = ORPHEUS `η` (radial); BMC/Hébert `η` = ORPHEUS `ξ` (azimuthal).
2. **Lathrop writes the weight as `Δμ_m`.** Read as a weight, not as a cell width
   in some other variable — although, by Eq. (C) itself, it *is* both. His level
   index and BMC's differ (`n` vs `p`) from Hébert's.
3. **The sign in front of the redistribution term is NOT universal**: BMC Eq. (49)
   carries `+[α₊ψ₊ − α₋ψ₋]/(r w)`, Hébert Eq. (3.400) carries `−(1/W)[α₊ψ₊ − α₋ψ₋]`.
   Hence BMC's `α_{m+1/2,n} = α_{m−1/2,n} − w μ` versus Hébert's
   `α_{p,q+1/2} = α_{p,q−1/2} + W μ`. **These are the same physics with opposite α
   sign conventions.** Check ORPHEUS's α sign against ORPHEUS's own streaming-term
   sign, not against a quoted recursion.

---

## 11. Zotero

Not consulted this pass — the brief scoped the task to the local folder ("ALL
LOCAL. Do NOT search online."), and all three sources were present with sidecars.
**No user annotations were checked**, so no annotation-based notation oracle is
available for this ruling. If a follow-up wants the user's own highlights on BMC
Eq. (42)/(74) or Hébert (3.398), that is a separate, cheap Tier-1 query.

---
---

# Follow-up: the transferable condition

Added 2026-08-11 in response to the internal-tension query. Everything below is
either scan-verified against the same three local PDFs or an explicitly marked
derivation. **One recommendation from §9 above is REFUTED here and corrected in
place — see §F0 and §F6.**

## F0. Corrections to the earlier sections (read first)

| earlier claim | status |
|---|---|
| §1 / §9(1): "(C)'s divergence is most likely the missing per-level weight renormalisation" | ⛔ **REFUTED** `[M]` F6. The renormalisation was applied correctly and (C) still puts τ outside `[0,1]`. The real cause is **node/partition incompatibility** — see F6. |
| §9(2): "β (BMC Eq. 75) is the adjudicator the MMS cannot be" | ⚠ **TRUE FOR THE SPHERE, FALSE FOR THIS CYLINDER RULE.** `[M]` F7: on a single polar level with equispaced-ω nodes and equal weights, β is **round-off for every convention including τ ≡ ½** — it is Mode-12 blind by the same oddness that Lathrop invokes in Eq. (57). A different diagnostic is needed; F7 gives it. |
| §4 item 2: "BMC's S₂ gives τ₁ ≈ 0.4226 < ½" | ✅ stands, and is much stronger than one value: `[M]` at **S₈ Gauss–Legendre, FOUR of the eight Morel–Montry τ are below ½** — `τ = [0.3923, 0.4591, 0.4809, 0.4942, 0.5058, 0.5191, 0.5409, 0.6077]`. A `max(0.5,·)` clamp corrupts **half** the values. |

## F1. THE TENSION RESOLVED — resolution (i), with the mechanism proved

**BMC's β and Lathrop's β are different objects that share a letter.** Side by side:

| | BMC (2010) Eq. **(41)/(75)**, p. 8 & 13 | Lathrop (2000) Eq. **(25)**, p. 8 |
|---|---|---|
| definition | `β_sph = ½ Σ_m η_m [α_{m+½}η_{m+½} − α_{m−½}η_{m−½}]` | `α_{m+½} = 1 − η_{m+½}² + β_{m+½}` |
| type | **one scalar** for the whole quadrature | **a sequence**, one per cell edge |
| what it measures | the coefficient of the O(ε) **contamination term in the second-order current** `J⁽²⁾ = −(1/3σ_t)∇φ⁽¹⁾ + (β/rσ_t²)∂φ⁽⁰⁾/∂r` | α's **pointwise additive departure from its exact geometric target** `1 − η_edge²` |
| vanishes when | the closure τ is the Morel–Montry one (Eq. 42/74) | `δ ≡ 0`, i.e. `τ ≡ ½` (Eq. 26) |
| scan status | `[VERIFIED-ON-SCAN p.8, p.13]` | `[VERIFIED-ON-SCAN p.8]` |

They vanish under **mutually exclusive** conditions. That is the whole of the
apparent contradiction: the coordinator read "β = 0" as one predicate when it is two.

### The mechanism — and the proof that it is the right reading

BMC's β is written with symbols `η_{m±½}` that look like the quadrature's
cumulative-weight edges. They are not. They are **the edge cosines that the
SCHEME's own τ implies**, through BMC Eq. (43) `[VERIFIED-ON-SCAN p.8]`, whose
role BMC state explicitly: Eq. (43) is what makes the flux closure Eq. (15)
*exact* when the angular flux is linear in the cosine. Solving Eq. (43) forward
gives the scheme's **effective edge march**

```
ν_{1/2} = −1 ,      ν_{m+1/2} = ( η_m − (1 − τ_m) ν_{m−1/2} ) / τ_m
```

and `β` is Eq. (41) evaluated on `ν`, with `α` still taken from the true
weights (Eq. 11).

`[M]` **This reading reproduces BMC's Table I to every printed digit** (sphere,
Gauss–Legendre; the table body was recovered from the OCR cache per L-008):

| N | step, printed | step, computed | diamond, printed | diamond, computed | Morel–Montry |
|---|---|---|---|---|---|
| 2 | 7.698004e-01 | **7.698004e-01** | 2.06E-01 | **2.06e-01** | 0.000e+00 |
| 4 | 4.037247e-01 | **4.037247e-01** | −3.57E-03 | **−3.57e-03** | 0.000e+00 |
| 8 | 2.164258e-01 | **2.164258e-01** | −4.57E-05 | **−4.57e-05** | −8.33e-17 |
| 10 | 1.755520e-01 | **1.755520e-01** | −1.21E-05 | **−1.21e-05** | 1.94e-16 |

At S₂ the picture is transparent: the quadrature's cumulative-weight edges are
`(−1, 0, +1)`; the Morel–Montry τ = `(0.42265, 0.57735)` reproduces them exactly,
so β = 0. The diamond's implied edges are `(−1, −0.15470, +1.30941)` — the interior
edge misses `μ_{3/2} = 0` (BMC's Eq. 46 requirement) and the march **overshoots the
physical range by 31 %**. That is the whole of β.

⟹ **Lathrop and BMC do not disagree.** Lathrop's §IV analysis of Eq. (18)
substitutes the *exact* linear flux at the *true* edges (`ψ_{m+1/2} = ψ(r, μ_{m+1/2})`),
and under that substitution his Eq. (57) shows the α-terms are sums of **odd**
functions against an even weight and vanish identically for any symmetric
quadrature — leaving `c = ⅔` as the only condition at that order. `[M]` My
reproduction confirms this: with the true cumulative-weight edges, β is round-off
for every N. The scheme-dependence enters only when the closure fails to realise
`ψ_{m+1/2} = ψ(r, μ_{m+1/2})`, which is exactly what Eq. (43) demands and only the
Morel–Montry τ delivers.

### Why Lathrop's δ=0 sentence is not about β at all

*"δ = 0 … would not have `c = ⅔` and hence would not give the correct diffusion
limit"* (`[VERIFIED-ON-SCAN p.12→13]`) is a statement about the **NODES that δ = 0
implies**: setting δ = 0 pointwise *defines* `μ_m = μ̄_m`, the cell midpoints, and
midpoint nodes fail the quadrature moment condition. `[M]` verified:

| N | Gauss–Legendre `c` | equal-interval MIDPOINT `c` | deviation |
|---|---|---|---|
| 2 | 0.666666666667 | 0.500000000000 | −1.667e−01 |
| 4 | 0.666666666667 | 0.625000000000 | −4.167e−02 |
| 8 | 0.666666666667 | 0.656250000000 | −1.042e−02 |
| 16 | 0.666666666667 | 0.664062500000 | −2.604e−03 |

(deviation `= −h²/6` with `h = 2/N`; it →0 but is never 0.)

⟹ **The two conditions live at two different ORDERS of the asymptotic expansion,
and this is BMC's own framing** (`[VERIFIED-ON-SCAN p.8]`: the leading-order flux is
accurate for *all* angular discretisations; Eq. (41) is what step and diamond fail,
"resulting in a first-order error in the diffusion limit"):

- **`c = ⅔` — a QUADRATURE condition, LEADING order.**
- **`β = 0` / Eq. (43) — a CLOSURE condition, FIRST order.**
- **Lathrop's `β_{m+½}` — neither; it is a bookkeeping identity for α.**

## F2. Q1 — the rule-agnostic diffusion-limit predicate

Lathrop derives `c` explicitly, which is exactly what was asked for.
`[VERIFIED-ON-SCAN p.12]` — Eq. (53) and Eq. (54): substituting a flux affine in
the marched cosine into the discrete system and taking the `Δμ_m`-weighted sum gives

```
d( r² · (3c/2) ψ_1 ) / (r² dr)  +  σ ψ_0  =  S_0 ,        c ≡ Σ_{m=1..N} η_m² Δη_m
```

and the P₁ pair (his Eqs. 52a/52b) is recovered **iff `3c/2 = 1`**, i.e.

> **(P1) `c ≡ Σ_m w_m η_m² = ⅔`** — Lathrop Eq. (29), which he names *the diffusion
> condition*. `[VERIFIED-ON-SCAN p.8, p.12]`

The second equation (his Eq. 59 = 52b) comes out right **unconditionally**, because
the α-terms reduce (his Eq. 57) to sums of odd functions. So (P1) is the entire
leading-order requirement. The remaining requirements, as a checkable predicate on
*any* 1-D angular march (nodes `η_m`, weights `w_m`, edges `e_{m±½}`, factors `τ_m`):

> **(P0) normalisation + odd closure.** `Σ_m w_m = |range|` and `Σ_m w_m η_m = 0`.
> The second is what makes the α recursion close (`α_½ = α_{M+½} = 0`) — stated by
> Lathrop `[VERIFIED-ON-SCAN p.12]` and by BMC Eq. (11)/(50).
>
> **(P1) `Σ_m w_m η_m² = ⅔ · Σ_m w_m / |range|`** (i.e. `⅔` under the `Σw = 2`
> normalisation). Leading-order diffusion limit. **A property of the QUADRATURE
> alone — no edges, no τ.**
>
> **(P2) `η_m = τ_m e_{m+½} + (1 − τ_m) e_{m−½}` for every m** — BMC Eq. (43) =
> Lathrop Eq. (23). First-order diffusion limit. **This is the condition on τ that
> is independent of how the edges were built** (see F5). Its scalar consequence is
> `β = 0`; its pointwise form is strictly stronger (F7).
>
> **(P3) `η_m ∈ [e_{m−½}, e_{m+½}]` for every m** — Lathrop's own gloss on Eq. (23):
> *"|δ| ≤ 1 as long as `μ_m` is in `Δμ`"* `[VERIFIED-ON-SCAN p.8]`. Equivalent to
> `τ_m ∈ [0,1]`. **Not an accuracy condition — a well-posedness one**: the closure
> is inverted as `ψ_{m+½} = ψ_m/τ_m − ((1−τ_m)/τ_m)ψ_{m−½}` (BMC Eq. 16 / 54), so
> `τ → 0` blows the recursion up and `τ < 0` flips its sign. **This is the predicate
> ORPHEUS's clamp was groping for, and clamping is the wrong repair** — a τ outside
> `[0,1]` means the *partition is incompatible with the nodes*, and the partition is
> what must change.
>
> **(P4) α satisfies the constant-flux recursion**: `α_{m+½} = α_{m−½} − k·w_m η_m`
> (`k = 2` sphere / `k = 1` cylinder-per-level), seeded from the geometric value at
> the range end. Hébert Eq. (3.397)/(3.398)/(3.424) `[VERIFIED-ON-SCAN p.72]`,
> BMC Eq. (11)/(50), Lathrop Eq. (19)/(58) `[VERIFIED-ON-SCAN p.12]`. A **first
> moment** condition — logically independent of (P2), which is a zeroth-order
> interpolation condition.

**Nothing in (P0)–(P4) mentions cumulative-weight edges.** BMC Eq. (12)/(52) is one
*solution* of the system, not the system. That is the answer to Q1 as posed.

## F3. Q2 — the cylinder's exactness target, and whether κ−1 is β

**The target: CONFIRMED, and it is Hébert's, verbatim.** `[VERIFIED-ON-SCAN p.72]`

Hébert defines `α_{p,q±1/2} ≡ W_p · η_{p,q±1/2}` under Eq. (3.399) — in ORPHEUS
symbols, **α = (polar weight) × (azimuthal/tangential cosine at the boundary)**,
which is exactly the form ORPHEUS's T3 theorem asserts. Its exactness condition is
Eq. (3.397)/(3.398), which in ORPHEUS symbols reads

```
w_m · η_m  =  w_gl · ( ξ_{m+1/2} − ξ_{m−1/2} )                    [Hébert (3.398)]
```

`[MY DERIVATION]` For equispaced ω with nodes at cell centres and edges at
`ω_m ∓ Δω/2`, using `ξ = sinθ sin ω`, `η = sinθ cos ω`:
`ξ_{m+½} − ξ_{m−½} = 2 sinθ · cos ω_m · sin(Δω/2)`, so the condition becomes
`w_m cos ω_m = w_gl · 2 cos ω_m sin(Δω/2)`. **The `cos ω_m` cancels identically**,
leaving `w_m = 2 w_gl sin(Δω/2)` for every m — which is precisely ORPHEUS's Q3
derivation, nodes untouched, uniform `1/κ` rescale. **Q3's algebra is confirmed.**

### Is `κ − 1` playing the role of Lathrop's β? — **YES, and that is why it is NOT the thing to drive to zero**

**CONFIRMED structurally.** Lathrop's `β_{m+½} := α_{m+½} − (1 − η_{m+½}²)` is
*α's additive departure from its exact geometric target*. ORPHEUS's α is `κ` times
its exact tangential target, so the additive departure is `(κ−1)·T_m`. Same object,
multiplicative vs additive bookkeeping. `[MY DERIVATION, structural]`

⛔ **But the rider is the whole point.** Lathrop Eq. (26) makes his β proportional
to δ, so **`β ≡ 0 ⟺ δ ≡ 0 ⟺ τ ≡ ½`** — the scheme with the *wrong* diffusion limit.
So confirming the analogy does **not** license driving `κ → 1`. Driving `κ → 1` is
driving the scheme toward angular diamond.

⚠ `[M]` **Arithmetic discrepancy in Lathrop Eq. (26), flagged not resolved.** From
his own Eq. (19) (`α_{m+½} − α_{m−½} = −2μ_m Δμ`, `[VERIFIED-ON-SCAN p.8]` incl. the
sentence taking `α_½ = 0 = 1 − μ_½²`) and Eq. (23), I get
`β_{m+½} − β_{m−½} = −2Δμ(μ_m − μ̄) = −δ Δμ²`, where the page prints `−δΔμ²/2`.
Factor 2. The **structure** — β linear in δ, zero iff δ ≡ 0 — is unaffected and is
what F1/F3 rest on. Do not transcribe the constant without re-deriving it.
(L-010 class: this is a *printed*-page reading, not an OCR slip.)

## F4. Q3 — do the sources confront the two-weights conflict? **No, because they never impose a geometric edge**

**Answer: the conflict cannot arise in their rule family, because the EDGE is the
free variable and the weight is given.** Both sources are explicit:

- **Hébert** takes `W_{p,q}` from the quadrature and *defines* `η_{p,q±½}` by the
  recursion (3.398). The half-angle symbols `ω_{p,q±½}` appear in Eq. (3.393) only
  as labels for where the half-angle flux lives; **no equation in §3.9.3 ever
  computes an ω half-angle.** `[VERIFIED-ON-SCAN p.72]`
- **BMC** likewise define the edges by a weight recursion (Eq. 12/52), and in their
  S₂ worked example they treat the interior edge as an **unknown determined by
  β = 0** (Eqs. 44→46 give `μ_{3/2} = 0`), then get τ from Eq. (42).
  `[VERIFIED-ON-SCAN p.8]`

⟹ **Neither source ever asks the edge to be a geometric bisection of the ω-cell.**
ORPHEUS's scheme does the opposite — it fixes the edges at the arc half-angles and
would then have to move the weights. The sources' option (**keep the nodes and
weights; let the edges float to whatever the recursion says**) preserves T1's
`Σw = π` trapezoid exactness *and* achieves α-exactness simultaneously, because
α-exactness is then a *definition* rather than a constraint.

`[MY DERIVATION]` And here is what that costs, in ORPHEUS's own quantity: run
Hébert (3.398) on ORPHEUS's nodes and trapezoid weights,
`ξ̃_{m+½} = ξ̃_{m−½} + Δω·sinθ·cos ω_m`, seeded `ξ̃_{−½} = 0`. Each increment is
`κ` times the geometric increment `2 sinθ cos ω_m sin(Δω/2)`, so

```
ξ̃_{m+1/2}  =  κ · ξ(ω_{m+1/2})          for every m,   and  ξ̃_end = κ·0 = 0  ✓
```

⟹ ⭐⭐ **ORPHEUS's production α ALREADY IS Hébert's α, exactly.** T3 is a
rediscovery of Eq. (3.398). `κ` is not an error inside α; it is the ratio between
**the recursion-defined edge — which is what α *is*, by definition, in both sources
— and the geometric arc-half-angle edge, which no source uses.** There is nothing
to fix in α, and the weights must not be rescaled: rescaling would break (P0)
(`Σw = |range|`, hence the collision term and T1) to buy exactness against a target
the literature does not adopt.

## F5. Q4 — the τ condition that does not care how the edges were built

**It exists, it is stated twice in the corpus, and it is the most transferable thing
in this whole document.**

> **BMC Eq. (43)** `[VERIFIED-ON-SCAN p.8]` ≡ **Lathrop Eq. (23)**
> `[VERIFIED-ON-SCAN p.8]`, with `τ = (1+δ)/2`:
> ```
> η_m  =  τ_m · e_{m+1/2}  +  (1 − τ_m) · e_{m−1/2}
> ```

In words, rule-agnostically: **τ_m is the barycentric coordinate of the ordinate's
own marched cosine between its two cell-boundary cosines.** Equivalently, and this
is the form to reason with: **τ_m is the unique weight for which the flux closure
`ψ_m = τ_m ψ_{m+½} + (1−τ_m) ψ_{m−½}` is EXACT whenever ψ is affine in the marched
cosine.** BMC say exactly this under Eq. (43), and Lathrop's Eq. (23)→(24) pairing
is the same statement written as a definition.

Three properties that make it the right target:

1. **It is a condition on the triple `(η_m, e_{m−½}, e_{m+½}, τ_m)` jointly**, and
   is silent about the provenance of the edges. Any partition may be used; τ is
   then *determined, never chosen*.
2. **τ is never a free parameter.** BMC: the τ of Eq. (42) are "uniquely defined by"
   Eq. (43). Lathrop's lineage note (`[VERIFIED-ON-SCAN p.12]`) is the same point
   historically: Morel–Montry keep Gauss–Legendre and **derive** δ; Reed–Lathrop
   instead choose the *cosines* and let the quadrature follow.
3. **Its scalar shadow is `β = 0`, and the shadow is strictly weaker.** `[M]`
   Perturbing `τ_1` at N = 8 and re-solving `τ_N` so the implied edge march still
   closes exactly at +1 leaves β non-zero (`ε = 1e−3 → β = −1.42e−06`;
   `5e−2 → −6.3e−05`). So neither `β = 0` nor edge-closure is a substitute for the
   pointwise condition. Use (P2) pointwise.

## F6. ⛔ Why (C) really diverges on ORPHEUS's rule — my §9 hypothesis REFUTED

`[M]` computed for a single polar level, equispaced ω, nodes at cell centres,
weights renormalised **correctly** to `Σ w̄ = 2 sinθ` (i.e. BMC Eq. 52 applied as
written). Values are `τ_m`, and `η/sinθ` edges:

| M | convention | τ range | node inside its own cell? |
|---|---|---|---|
| 8 | **(C)** cumulative-weight (uniform-η cells) | **[−0.3259, 1.3259]** | **NO** — fails at m = 2,3,6,7 |
| 16 | **(C)** | **[−1.1841, 2.1841]** | **NO** — worse |
| 32 | **(C)** | **[−2.8552, 3.8552]** | **NO** — worse still |
| 8 | (B) arc half-angle | [0.2524, 0.7476] | yes, all M |
| 16 | (B) | [0.2506, 0.7494] | yes |
| 32 | (B) | [0.2502, 0.7498] | yes |
| 8 | (A) chord midpoint | [0.2047, 0.7953] | yes |
| 32 | (A) | [0.2003, 0.7997] | yes |

⟹ **The renormalisation was never the problem.** (C) violates **(P3)**: with
equispaced-ω nodes, the ordinate's radial cosine `η_m = sinθ cos ω_m` simply does
not lie inside the m-th *uniform-η* cell. And the violation **grows with
refinement** — which matches the reported ladder exactly (finite at `n_phi = 8`,
`nan` at 16 and 32).

**The structural reason, and it is the concept that transfers.** `[MY DERIVATION]`
The η-measure of arc cell m is `2 sinθ · sin(ω_m) · sin(Δω/2)` — **proportional to
`sin ω_m`, not constant** — while the trapezoid weight *is* constant. Measured
ratio (measure / weight) at M = 8: `[0.30, 0.87, 1.30, 1.53, 1.53, 1.30, 0.87, 0.30]`.

> ⭐ **BMC Eq. (52) is not a law. It is the statement that in THEIR quadrature
> family the ordinate weight happens to equal the cell's η-measure.** For an
> equispaced-ω rule that identity is false by up to 5×, so Eq. (52) is simply not
> available — and forcing it produces a partition whose cells do not contain their
> own nodes. **(C) is the wrong family member for ORPHEUS's rule, and no
> renormalisation repairs that.**

### And (A) is (B) in disguise

`[M]` to 1e−16 at M = 8, 16, 32: **ORPHEUS's (A) interior edges are exactly
`cos(Δω/2) × (B)`'s arc edges**, with only the two end cells stretched to reach
`∓sinθ`. So (A) and (B) are the *same* family — (A) is (B) contracted by
`cos(Δω/2) = 1 − Δω²/8 + …` with inconsistent end cells. That is why they measured
so similarly, and it makes **(B) the clean member**: it is the same partition
without the O(Δω²) contraction and without the endpoint inconsistency.

`[MY DERIVATION]` (B) has a closed form, verified to 6.7e−16 at M = 8/16/32:
```
τ_m  =  ½  +  ½ · cot(ω_m) · tan(Δω/4)
```
— symmetric about `ω = π/2` (where τ = ½ exactly), departing only near the march
ends, and → ½ as `Δω → 0` at fixed ω. Note this satisfies (P3) for all M
(`max|τ − ½| → ¼` from below) and reduces to the diamond in the interior, which is
the qualitative behaviour Lathrop reports for Gauss–Legendre δ
(`[VERIFIED-ON-SCAN p.8]`: δ small except near the ends of the range, decreasing
with order).

## F7. ⚠ The oracle: β is BLIND on this rule — use (P2)/(P3) pointwise instead

`[M]` BMC Eq. (75) evaluated per polar level on the equispaced-ω march
(α from Eq. 50 with ORPHEUS's weights, edges from each scheme's implied march):

| M | (B) arc | (C) cumul-weight | (A) chord | **τ ≡ ½** | α closes? |
|---|---|---|---|---|---|
| 8 | −1.7e−16 | −8.3e−17 | −1.7e−16 | **−2.8e−16** | yes, 1.7e−16 |
| 16 | −2.2e−16 | 4.4e−16 | −1.1e−16 | **−2.2e−16** | yes |
| 32 | −8.3e−17 | 5.5e−15 | −8.3e−17 | **−8.3e−17** | yes |

**β is round-off for every convention, including τ ≡ ½.** The cause is the same
oddness Lathrop invokes in his Eq. (57): an equispaced-ω march with equal weights
is symmetric about `ω = π/2`, `η → −η` maps the node set to itself, and every term
cancels pairwise. So on ORPHEUS's rule **β cannot see the closure choice at all** —
it is a Mode-12 blind functional here, exactly as the MMS is. (It is *not* blind on
the sphere with Gauss–Legendre, where it reproduces Table I to all digits; the
blindness is a property of this node/weight symmetry, not of β.)

⟹ **Use these instead, in this order — all cheap, none needing a solve:**

1. **(P3) `τ_m ∈ [0,1]` for all m, at every `n_phi`.** Kills (C) outright and is
   the well-posedness gate the clamp was standing in for. Assert the *predicate*,
   never clamp the value.
2. **Edge-march closure**: the τ-implied march `ν` (F1) must land exactly on
   `+sinθ`. `[M]` (A), (B), (C) all close exactly; **τ ≡ ½ does not** —
   `ν_end/sinθ = 1.0392 (M=8) → 1.0097 (16) → 1.0024 (32)`, an O(Δω²) overshoot of
   the physical range. **This is the one diagnostic that separates a derived τ from
   the diamond on ORPHEUS's own rule**, and it is one line of arithmetic.
3. **(P2) pointwise residual**: `max_m |η_m − τ_m e_{m+½} − (1−τ_m)e_{m−½}|`. Zero
   by construction for any derived τ; O(Δω²) for τ ≡ ½. This is the *definition*,
   so it is the honest gate.
4. **(P1) `Σ w η² = ⅔ · Σw/|range|`** on the polar rule — a pure quadrature check,
   independent of everything above.

## F8. Bottom line for the campaign

- **Keep the equispaced-ω folded arc and the trapezoid weights.** Do not rescale
  the weights (F4): that would break (P0)/T1 to chase a target no source adopts,
  and ORPHEUS's α already equals Hébert's α exactly.
- **Use the arc partition (B)** — not on BMC's authority (they prescribe (C) for
  *their* nodes) but because (B) is the partition that satisfies **(P3)** for
  ORPHEUS's node set, and (P2) then *determines* τ with the closed form in F6.
  (A) is (B) with an O(Δω²) contraction and broken end cells; prefer (B).
- **Delete the clamp; assert (P3) as a predicate.** τ below ½ is normal — at S₈ on
  the sphere, half the Morel–Montry values are.
- **Do not use β as the gate on this rule** (F7). Use closure + the (P2) residual.
- ⚠ **What no source settles for ORPHEUS's rule**: whether some *third* partition
  (edges floated from a second exactness condition, as Hébert floats the azimuthal
  ones) would beat (B). The corpus has no cylinder analogue of BMC's Eq. (52)
  derivation for non-cosine-distributed nodes. That is an original-work gap, not a
  reading gap — and **Morel & Montry (1984) TTSP 13(5):615 is still the one
  unacquired primary source** where a general statement of the closure condition
  could live.
