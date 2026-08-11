# q64 attempt 2 — published check of two 1-D cylindrical S_N claims

Deliverable of a literature pull, 2026-08-11. Local library only unless noted.
Written INCREMENTALLY (content-filter discipline). Verbatim quotes kept SHORT.

**Notation of the BRIEF (= ORPHEUS):** ω azimuthal angle from the local outward
radial unit vector `e_r`; θ polar from `z`; `η = Ω·e_r = sinθ cosω` (radial cosine);
`ξ = Ω·e_φ = sinθ sinω` (azimuthal cosine); `μ_z = Ω·e_z = cosθ`.

⚠⚠ **THE MASTER NOTATION CLASH — read before every finding below.**
Hébert (and BMC, and Lathrop) use the SAME three letters for a DIFFERENT
assignment. Translation table, valid for Hébert 2009 Ch. 3 and BMC 2010:

| role | ORPHEUS / brief | Hébert 2009 | BMC 2010 (R-Z) | Lathrop 2000 (sphere) |
|---|---|---|---|---|
| **radial** cosine `Ω·e_r` | `η` | **`μ`** | **`μ`** | `μ` (is the only cosine) |
| **azimuthal** cosine `Ω·e_φ` | `ξ` | **`η`** | **`η`** | — |
| **axial/polar** cosine `Ω·e_z` | `μ_z` | **`ξ`** | **`ξ`** | — |
| azimuthal angle from `e_r` | `ω` | `ω` | `ω` (implicit) | — |

⟹ **Hébert's `μ` is the brief's `η`, and Hébert's `η` is the brief's `ξ`.** Every
equation below is quoted in the SOURCE's letters and then translated. A reader who
skips this table will read the cylinder's redistribution coefficient as the radial
cosine, which is the exact opposite of what it is.

---

# CLAIM 1 — the conservative streaming form in 1-D cylindrical (r)

> `Ω·∇ψ = (η/r) ∂(rψ)/∂r − (1/r) ∂(ξψ)/∂ω`

## 1.1 — Hébert 2009, Ch. 3, printed p. 92, **Eq. (3.157)** — `SUPPORTS` (exact, verbatim)

`[VERIFIED-ON-SCAN PDF p. 26 = printed p. 92]`

The printed equation, in HIS letters (`μ`=radial, `η`=azimuthal, `ξ`=polar):

```
(μ/ρ) ∂/∂ρ[ρ φ(ρ,ξ,ω)]  −  (1/ρ) ∂/∂ω[η φ(ρ,ξ,ω)]  +  Σ(ρ) φ(ρ,ξ,ω)  =  Q
```

Introduced by one sentence: *"The conservative form of Eq. (3.155) is written"*.

Translated into the brief's letters this is **character-for-character the claim**:
`(η/r) ∂(rψ)/∂r − (1/r) ∂(ξψ)/∂ω`. Same `1/r` prefactor on both terms, same minus
sign on the azimuthal term, no `sinθ` anywhere, no `1/sinθ`.

**VERDICT: `SUPPORTS`.** This is the equation, published, in the local library.

## 1.2 — Hébert 2009, printed p. 91, **Eq. (3.155)** — the non-conservative twin `SUPPORTS`

`[VERIFIED-ON-SCAN PDF p. 25 = printed p. 91]`

```
μ ∂φ/∂ρ  −  (η/ρ) ∂φ/∂ω  +  Σ(ρ) φ  =  Q
```

⟹ brief's letters: `η ∂ψ/∂r − (ξ/r) ∂ψ/∂ω`. **Exactly the alternative form the
brief offers.** Both forms are printed, two lines apart, and (3.157) is explicitly
derived from (3.155).

**Consistency check (L-010 (b), non-vacuity) — the two forms are mutually
consistent, and Hébert states the identity that makes them so.** Same page, one
line under Eq. (3.152): *"Eqs. (3.152) can also be used to show that `∂η/∂ω = μ`"*.
So expanding (3.157): `(μ/ρ)(ρ∂φ/∂ρ + φ) − (1/ρ)(μφ + η ∂φ/∂ω)` = `μ∂φ/∂ρ −
(η/ρ)∂φ/∂ω` — the two spurious `μφ/ρ` terms cancel **identically**. In the brief's
letters the identity is `∂ξ/∂ω = η`. ⟹ the conservative form is not an
approximation or a rearrangement with a leftover; it is exact.

## 1.3 — Provenance of the general 3-D form: Hébert Eqs. (3.147)–(3.151) — `SUPPORTS`

`[VERIFIED-ON-SCAN PDF p. 24–25]`. 3-D cylindrical, before the 1-D reduction:

- (3.148): `dρ/ds = μ`, `dθ/ds = η/ρ`, `dz/ds = ξ`.
- (3.149): `dω/ds = −(1/ρ)√(1−ξ²) sinω = −η/ρ`.
- (3.151): `Ω·∇φ = (μ/ρ)∂(ρφ)/∂ρ + (η/ρ)∂φ/∂θ + ξ ∂φ/∂z − (1/ρ)∂(ηφ)/∂ω`.

The 1-D form drops `∂/∂θ` and `∂/∂z` by symmetry (3.153) and yields (3.155)/(3.157)
with **nothing else changed**. ⟹ the claimed 1-D form is a strict specialisation of
a published 3-D form, not a 1-D-only construct.

## 1.4 — Claim 1(b): the azimuthal FACE coefficient. Hébert **Eq. (3.393)** — `SUPPORTS`, with one caveat that is the whole of Claim 2

`[VERIFIED-ON-SCAN — pending page read, see §1.4v]` Sidecar text of (3.393):

```
∫_{Ω_n} d²Ω ∂/∂ω[η φ]  =  W_p [ η_{p,q+1/2} φ(ρ, ξ_p, ω_{p,q+1/2})
                                − η_{p,q−1/2} φ(ρ, ξ_p, ω_{p,q−1/2}) ]
```

- The coefficient multiplying the half-angle (face) flux is **`η` at that face and
  nothing else** — i.e. the brief's `ξ` at the face. **No extra geometric factor.**
- The one factor that IS present, `W_p`, is the **polar-level weight** — the measure
  of the `ξ`-integration for level `p`, not a geometry factor. It divides out: the
  balance equation (3.400) carries the face term as
  `−(1/W_{p,q})[α_{p,q+1/2} φ_{p,q+1/2} − α_{p,q−1/2} φ_{p,q−1/2}]` with
  `α_{p,q±1/2} ≡ W_p η_{p,q±1/2}` (3.399).
- ⚠ **CAVEAT, and it is the crux of Claim 2:** Hébert does **not** evaluate
  `η_{p,q±1/2}` geometrically as `√(1−ξ_p²) sin ω_{p,q±1/2}`. He *defines* it by the
  recursion (3.398) `η_{p,q+1/2} = η_{p,q−1/2} + (W_{p,q}/W_p) μ_{p,q}`, seeded by
  (3.396) `η = 0` at the level's end point. So the FORM of Claim 1(b) is exactly as
  the brief states it; the VALUE of `ξ_face` in a discretised scheme is a
  recursion output, and only coincides with the geometric arc value if the
  quadrature weight equals the cell's `η`-measure.
- (3.158), printed p. 92: `∫_{4π} d²Ω ∂/∂ω[ηφ] = 0` — the redistribution term is a
  pure telescoping face flux with zero total. This is the continuous statement whose
  discrete image IS the recursion (3.397)/(3.398).

## 1.5 — Claim 1(c): the ω sign/ordering convention — `SUPPORTS` ORPHEUS's storage

Three independent statements in Hébert, all `[VERIFIED-ON-SCAN]`:

1. **ω is measured from `e_ρ`** — Fig. 3.5 + Fig. 3.6 (printed p. 90–91) show ω as
   the angle between the projected direction and `e_ρ`, and (3.152) `μ = √(1−ξ²)cosω`
   makes it unambiguous: **ω = 0 ⟹ maximally outward; ω = π ⟹ maximally inward.**
   Same convention as the brief.
2. **ω DECREASES along the particle path.** Eq. (3.149) is negative, and Hébert says
   why in one clause: *"the negative sign indicates that the azimuthal angle ω
   decreases as the particle goes forward"* (printed p. 90).
3. **The level is traversed from ω = π to ω = 0.** The sweep description (printed
   p. ≈107, sidecar line 2921 + Fig. 3.27) starts the level at ω = π — the
   `μ < 0` (inward) end — and Eq. (3.407) is introduced with *"where ω = π since the
   particle is moving toward the central axis"*. The march then proceeds through
   decreasing ω, i.e. **increasing radial cosine μ**.

⟹ **ORPHEUS storing the level η-ASCENDING (ω-DESCENDING, π → 0) is Hébert's own
traversal order.** `SUPPORTS`.

⚠⚠ **BUT SEPARATE "STORAGE INDEX ORDER" FROM "SWEEP ORDER" — Hébert's two are
OPPOSITE, and BMC's two agree.** This is the single most confusable point in Claim 1(c):

| | storage index runs | sweep runs |
|---|---|---|
| **Hébert 2009** | `q` **increases with ω** (Eq. 3.370, §2.6) ⟹ ω: 0 → π | **ω: π → 0** (3.407 + Fig. 3.27) — i.e. **against** the index |
| **BMC 2010** | `m` **increases with the radial cosine μ** (Eq. 52, seed `−√(1−ξ²)`) ⟹ ω: π → 0 | ω: π → 0 — **with** the index |
| **Carlson–Lathrop 1965** | `μ_m` increases on `[−1,1]` (printed p. 12) ⟹ ω: π → 0 | with the index |
| **ORPHEUS** | η-ascending ⟹ ω: π → 0 | ω: π → 0 |

⟹ **ORPHEUS matches BMC and Carlson–Lathrop on BOTH counts.** It matches Hébert's
*sweep* but is the *reverse of Hébert's index*. Consequence: any formula lifted from
Hébert's `q`-indexed equations must have its `q±½` labels swapped — and §2.7 shows this
flips the sign of a term in the closed-form τ.

## 1.6 — Pomraning 1989 NSE 101(4):330–340 — `SILENT` (and explicitly so)

`[sidecar]` The paper derives a general-geometry streaming operator (its Eqs. (68),
(75) carry the angular-derivative coefficients as reciprocal radii of curvature) but
**declines the standard cylindrical case**: printed §V says the usual cylindrical
treatment *"is distinctly different from our considerations in this geometry"* and
that a geometry-tailored treatment "will probably be simpler and more natural",
because a finite cylinder's surface is convex but **not smooth** (the end-cap
corners make θ, φ non-smooth on the boundary).

**VERDICT: `SILENT` on the claim.** Do NOT cite Pomraning 1989 for the 1-D
cylindrical streaming form. What it IS good for: it names the provenance chain —
its refs. 8–10 are "the archival literature" for the standard cylindrical/spherical
streaming operator, and it states that *"there seems to be no easily accessible
derivation of these results"*, the only one found being a LANL report by **Lee**
(its ref. 11). That is a citable statement that Hébert's derivation is unusually
explicit for this operator.

## 1.7 — ⭐ Bell & Glasstone 1970, printed **p. 58** (display) and **p. 59 (Table 1.2)** — `SUPPORTS`, INDEPENDENT of Hébert

`[VERIFIED-ON-SCAN PDF pp. 76–77 = printed pp. 58–59]`

⚠ **NOTATION CLASH — B&G's `μ` is the AXIAL cosine, and their cylindrical azimuthal
angle is `χ`, not `ω`.** Printed p. 57: `μ = Ω̂·ẑ`; `χ` is the dihedral angle
locating `Ω` about `ẑ` relative to `r̂` (Fig. 1.16). They reserve `ω` for the
*spherical* azimuth. ⟹ `√(1−μ²) cos χ` = brief's `η`; `√(1−μ²) sin χ` = brief's `ξ`;
their `μ` = brief's `μ_z`. **A reader who maps B&G's `μ` onto the radial cosine gets
every coefficient wrong.**

**(a) The general 3-D conservation form** (printed p. 58, unnumbered display,
introduced by *"The following expressions give `Ω·∇Φ` in conservation form"*):

```
(√(1−μ²) cosχ / r) ∂(rΦ)/∂r  +  (√(1−μ²) sinχ / r) ∂Φ/∂φ
                              −  (1/r) ∂(Φ √(1−μ²) sinχ)/∂χ  +  μ ∂Φ/∂z
```

Drop `∂/∂φ` (axial symmetry) and `∂/∂z` (infinite cylinder) ⟹ in the brief's letters
**`(η/r)∂(rΦ)/∂r − (1/r)∂(ξΦ)/∂ω`**. Identical to Claim 1, term for term, sign for
sign, with **no `sinθ` factor and no `1/sinθ`**.

**(b) Table 1.2 row "Cylindrical (infinite cylinder, axial symmetry) `Φ(r, μ, χ)`"**
— the *non-conservative* twin, printed for exactly the 1-D case the brief asks about:

```
√(1−μ²) cosχ ∂Φ/∂r  −  (√(1−μ²)/r) sinχ ∂Φ/∂χ        [∫dΩ = ∫_{−1}^{1}dμ ∫_0^{2π}dχ]
```

⟹ brief's letters: **`η ∂Φ/∂r − (ξ/r) ∂Φ/∂ω`**. `SUPPORTS`, exactly.

**(c) B&G also state the CONSERVATION property that makes 1(b) true** (printed p. 59,
under Table 1.2): integrating the spherical conservation form over all directions
*"removes the last two terms, whereas the first three terms represent the components
of ∇·J"*. Same statement as Hébert (3.158), for the same reason. ⟹ the angular term
is a pure telescoping face flux with zero angular total.

**(d) The SIGN ASYMMETRY between the geometries is PHYSICAL, not conventional.**
B&G's cylindrical angular term is `−(1/r)∂(ξΦ)/∂χ`; their spherical one is
`+(1/r)∂[(1−μ²)Φ]/∂μ` (also Hébert (3.162), also Carlson–Lathrop (3-16)). The
cylinder is negative because the marched variable `ω` **decreases** along the path
(Hébert 3.149) while the sphere's `μ` **increases**. ⟹ never port a sign from the
sphere's redistribution term to the cylinder's.

## 1.8 — Carlson & Lathrop 1965, LA-3251-MS §3.1.4, printed pp. 11–14 — `SUPPORTS` the DISCRETE form of 1(b); `SILENT` on the closed-form cylinder PDE

`[sidecar; equations 3-9 … 3-18]`

The original S_N derivation, and *structurally* the source of 1(b):

- **Eq. (3-9)** — general one-way-curved conservative difference form, whose angular
  term is `(α_{m+½}N_{m+½} − α_{m−½}N_{m−½})/w_m`. Their words for what α is:
  *"the terms in α are introduced as flow terms from the 'edges' or 'surfaces' of the
  directional cell"* — the face flux `α_{m−½}N_{m−½}Δt` being *"directly analogous to
  flow through the geometric surface of the cell"*.
- **Eq. (3-10)** — the α recursion, from uniform flow:
  `α_{m+½} − α_{m−½} = − w_m μ_m (A_{i+1} − A_i)`.
- **Eqs. (3-11)/(3-12)** — conservation forces `α_{½} = α_{M+½} = 0`, hence
  `Σ_m w_m μ_m = 0` as *"an additional constraint on the angular quadrature
  coefficients"*.
- **Eq. (3-18)** — for a geometry *"curved in two directions"*, a second family
  `β_{m+½} − β_{m−½} = −w_m η_m (B_{j+1} − B_j)`, `β_{½}=β_{M+½}=0`.
  ⚠ **`β` here is a THIRD object wearing that letter** — neither BMC's scalar `β`
  (Eq. 41/75) nor Lathrop-2000's sequence `β` (Eq. 25). It is the *y*-curvature
  redistribution coefficient. Do not conflate.
- ⚠ **NOTATION: LA-3251 assigns `μ`↔x, `η`↔y, `ξ`↔z (Eq. 3-8), and a footnote
  (printed p. 17) says the letters PERMUTE for SPHERICAL geometry**
  (`ξ_m→μ_m, μ_m→η_m, η_m→ξ_m`). ⟹ for cylindrical, `μ` = radial, `ξ` = axial,
  `η` = azimuthal — the **same** assignment as Hébert and BMC, and the *opposite*
  of B&G.
- **`SILENT`** on the closed-form cylindrical PDE: the small-interval limit is taken
  only for the **sphere**, giving (3-16) `(μ/r²)∂(r²N)/∂r + (1/r)∂[(1−μ²)N]/∂μ + σN
  = S`. No cylinder limit is printed.
- ⭐ **Their ordering statement is ORPHEUS's**: the α's are non-negative and the flows
  interpretable *"with positive `w_m`, `A_{i+1} > A_i`, and direction sets ordered so
  that the `μ_m` increases uniformly on the interval [−1,1]"*. In cylindrical, `μ` is
  the radial cosine ⟹ **the level is traversed radial-cosine-ASCENDING.**
  Independent confirmation of Claim 1(c).

## 1.9 — ⛔ Stacey 2007 Ch. 9, printed p. 344, **Eq. (9.150)** — `SUPPORTS` in intent, but the PRINTED coefficient is WRONG (L-010 class)

`[VERIFIED-ON-SCAN PDF p. 41 = printed p. 344 — the OCR is FAITHFUL; the PRINT is the error]`

Stacey's notation (Eq. 9.149): `μ = cosθ = Ω·n_z` (axial, like B&G), and `φ` is *"the
angle in the x–y plane between the x–y projection of Ω and the radius vector r"* —
i.e. **exactly the brief's `ω`**. So `sinθ cosφ` = brief `η`, `sinθ sinφ` = brief `ξ`.

The equation as **printed**:

```
sinθ [ cosφ ∂ψ/∂r  −  (sinθ/4) ∂ψ/∂φ ]  +  Σ_t ψ  =  ...
```

⛔ **A typographical error in the published textbook**, provable two ways with no
external source:
1. **Dimensions.** `∂ψ/∂r` carries `1/length`; `∂ψ/∂φ` is dimensionless, and `sinθ/4`
   is dimensionless ⟹ the bracket adds `1/length` to a pure number. The printed
   equation is dimensionally inconsistent, hence cannot be right.
2. **The twin.** Factor `sinθ` out of B&G's Table 1.2 infinite-cylinder row:
   `sinθ[cosχ ∂Φ/∂r − (sinχ/r)∂Φ/∂χ]`. The intended Stacey coefficient is therefore
   **`sinφ/r`** — i.e. `4 → r` AND `sinθ → sinφ`, two slips inside one bracket.

⟹ **VERDICT: `SUPPORTS` once corrected**, but **this is the one local source that,
read literally, "writes the angular derivative term with a different coefficient
involving sinθ separately"** — exactly the failure mode the brief asked about. Do
**not** cite Stacey (9.150) as printed. If a theory page wants a second citation
beside Hébert, use **Bell & Glasstone Table 1.2**, which is correct as printed.

## 1.10 — `SILENT` sources checked and excluded (reason given so nobody re-checks)

| source | why `SILENT` |
|---|---|
| Ligou 1982 Ch. 8 | Integral-transport / collision-probability + `P_1` boundary theory. No `S_N` streaming operator in any geometry. |
| Pomraning 1989 | §1.6 — general geometry only, and *declines* the cylinder by name. |
| LA-3186, LA-4058 | Quadrature-set reports (level-symmetric ordinates/weights). No streaming operator. |
| Lathrop 2000 | 1-D **spherical** throughout, by title and body. No cylinder anywhere. |
| Sanchez–Ganapol 1983 | Integral-transport / singular-eigenfunction cylinder benchmark; no differential form with an ω-derivative. |

## 1.11 — CLAIM 1 SUMMARY

| sub-claim | verdict | best citation |
|---|---|---|
| `Ω·∇ψ = (η/r)∂(rψ)/∂r − (1/r)∂(ξψ)/∂ω` | **`SUPPORTS`** | Hébert 2009 **Eq. (3.157)** p. 92; B&G 1970 p. 58 display |
| non-conservative `η ∂ψ/∂r − (ξ/r)∂ψ/∂ω` | **`SUPPORTS`** | Hébert **Eq. (3.155)** p. 91; B&G **Table 1.2** infinite-cylinder row p. 59 |
| (b) azimuthal face coefficient = `ξ` at the face, no extra geometric factor | **`SUPPORTS`** in form; ⚠ the *value* is recursion-defined, not geometric | Hébert **(3.393)** + (3.399); Carlson–Lathrop **(3-9)/(3-10)** |
| (c) ω from `e_r`, DECREASES along the path, level traversed π → 0 | **`SUPPORTS`** ORPHEUS's η-ascending storage | Hébert **(3.149)** + **(3.407)** + Fig. 3.27; Carlson–Lathrop printed p. 12 |
| any source with a different angular coefficient | **only Stacey (9.150) — a print typo** | §1.9 |

## 1.12 — ⭐ BONUS: BMC 2010 **Eq. (48)** is a THIRD independent statement of Claim 1

`[VERIFIED-ON-SCAN PDF p. 9 = printed p. 156]`. Introduced as *"The R-Z transport
equation … is written in conservation form¹ as"* (their ref. 1 — a textbook — is
credited for the form, not derived here):

```
(μ/r) ∂/∂r rψ(r,z,ω,ξ)  −  (1/r) ∂/∂ω ηψ(r,z,ω,ξ)  +  ξ ∂ψ/∂z  +  σ_t ψ  =  (1/4π)σ_s φ + (1/4π)Q
```

BMC's letters are Hébert's (`μ`=radial, `η`=azimuthal, `ξ`=axial), so this is again
**exactly Claim 1** plus the axial term. Drop `ξ ∂/∂z` for the infinite cylinder ⟹
identical to Hébert (3.157). `SUPPORTS`.

⛔ **BUT THE SENTENCE DEFINING THE SYMBOLS IS A PUBLISHED TYPO — flag it, because it
inverts the very quantity Claim 1 is about.** Printed p. 156, immediately under
Eq. (48): *"where ξ = cos θ, μ = sin θ cos ω, η = sin θ cos ω in Eq. (48)"* —
**`cos ω` appears twice.** `[VERIFIED-ON-SCAN: the OCR is faithful; the PRINT is
wrong.]` The correct definition is **`η = sin θ sin ω`**, provable four ways with no
external appeal:
1. `μ` and `η` cannot both equal `sinθ cos ω` (they are two different components).
2. `μ² + η² + ξ² = 1` fails for the printed pair.
3. The conservative form requires `∂η/∂ω = μ` (Hébert's stated identity), true iff
   `η = sinθ sin ω`.
4. Their Fig. 2 draws `ω` as the azimuth of `Ω` from `e_r` toward `e_γ`.
⟹ A reader taking the printed definition literally would conclude the azimuthal-face
coefficient is the RADIAL cosine. **This is the single most dangerous notation trap
in the corpus for this claim.**

⚠ Two further published slips in the same section, all the same self-referential
family (see also §2.7): **Eq. (50)** prints `α_{m+1/2,n} = α_{m+1/2,n} − μ_{m,n}w_{m,n}`
(RHS must be `α_{m−1/2,n}`), and **Eq. (52)** prints
`μ_{m+1/2,n} = μ_{m+1/2,n} + w̄_{m,n}` (RHS must be `μ_{m−1/2,n}`).
`[VERIFIED-ON-SCAN pp. 9–10]`

⚠ **α SIGN.** BMC (50) is `α_{m+½} = α_{m−½} − μ_m w_m`; Carlson–Lathrop (3-10) is the
same sign; **Hébert (3.399) is `+`**. The difference tracks the sign in front of the
redistribution term in each author's balance equation (BMC/C-L: `+[…]/(rw)`;
Hébert: `−(1/W)[…]`). Check ORPHEUS's α sign against ORPHEUS's own streaming sign,
never against a quoted recursion.

---

# CLAIM 2 — how the sources derive τ for the CYLINDER, and in WHICH variable

## 2.1 — (a) THE cylinder formula exists: BMC 2010 **Eq. (74)**, printed p. 160 — barycentric in the **RADIAL COSINE**, not in ω. `SUPPORTS` ORPHEUS

`[VERIFIED-ON-SCAN PDF p. 13 = printed p. 160]`

```
τ_{m,n}  =  (μ_{m,n} − μ_{m−1/2,n}) / (μ_{m+1/2,n} − μ_{m−1/2,n})              (74)
```

introduced by: *"This term will be zero if we use the Morel and Montry weighting
factors in Eq. (53) analogous to those defined by Eq. (42)"*.

- BMC's `μ` **is the radial cosine** (`μ = sinθ cos ω`, printed under Eq. 48) ⟹ in the
  brief's letters, `τ_{m,n} = (η_m − η_{m−1/2}) / (η_{m+1/2} − η_{m−1/2})`, all three
  cosines taken **on one polar level `n`**.
- **⟹ the barycentric coordinate is in the RADIAL COSINE η. It is NOT in the azimuthal
  angle ω.** No source in the corpus takes it in ω. `SUPPORTS` ORPHEUS's
  implementation exactly, and **REFUTES** an ω-barycentric reading.
- The cells that define the edges are **radial-cosine cells on a level**: Eq. (52)
  `μ_{m+1/2,n} = μ_{m−1/2,n} + w̄_{m,n}`, `μ_{1/2,n} = −√(1−ξ_n²)`,
  `μ_{M+1/2,n} = +√(1−ξ_n²)`, *"where the weights in Eq. (52) are normalized on each
  level to sum to `2√(1−ξ_n²)`"*. `[VERIFIED-ON-SCAN p. 10]`
- ⚠ **Eq. (74) is verbatim identical in form to the sphere's Eq. (42)** — the paper
  says "analogous to". So ORPHEUS's code comment "BMC Eq. 43 = Lathrop Eq. 23" is
  citing the *sphere* equation for a *cylinder* implementation. **The cylinder has its
  OWN equation number and it should be cited: Eq. (74) (edges: Eq. 52; closure:
  Eq. 53).** Not an error of substance — an error of provenance.

**`SILENT`/absent elsewhere:** Hébert 2009 has **no τ at all** for either geometry —
his cylinder angular closure is the plain diamond, (3.406), and his sphere (3.431).
Lathrop 2000 and LA-3251/LA-3186/LA-4058 are sphere-only or quadrature-only.
⟹ **BMC Eq. (74) is the LONE cylinder τ statement in the local corpus.**

## 2.2 — (b) What exactness is claimed? **TWO different properties, both printed, and they are linked.** `SUPPORTS` both readings the brief offered

**(i) "Exact for a flux affine in the RADIAL cosine."** BMC printed p. 155, the
sentence immediately after Eq. (43):

> *"Equation (43) implies that Eq. (15) will exactly relate the cell-edge and
> cell-center fluxes when the angular flux assumes the linear form defined by
> Eq. (1)."*

with Eq. (15) = the weighted-diamond closure `ψ_m = τ_m ψ_{m+½} + (1−τ_m) ψ_{m−½}`,
Eq. (43) = `μ_m = τ_m μ_{m+½} + (1−τ_m) μ_{m−½}`, and **Eq. (1)** (printed p. 150,
their step 2 of the Galerkin procedure, introduced *"Assume a discrete angular flux
solution of the following linear form"*):

```
ψ_m  =  (1/4π) φ  +  (3/4π) J_r μ_m                                            (1)
```

where they gloss `J_r` as *"radial component of the current"*.

**(ii) "First-order asymptotic diffusion-limit consistency" = "preserves the Galerkin
diffusion approximation."** BMC printed p. 160, right after Eq. (74):

> *"Thus, if we use the weighting factors given in Eq. (74), we will achieve
> first-order consistency in the diffusion limit"*

and under Eq. (75): *"This is the same factor that Morel and Montry forced to be zero
to preserve the Galerkin diffusion approximation."* Plus the uniqueness claim
(printed p. 161): the MM scheme *"is the only method in the family of general weighted
diamond methods (which includes the step and diamond schemes) that forces the β factor
to be zero for any standard quadrature set."*

⭐ **The two are the same property seen from two sides, and BMC's own Eqs. (35)/(61)
say why**: the diffusion-limit first-order angular flux IS affine in the radial cosine
(sphere Eq. (35): `ψ_m^{(1)} = ½φ^{(1)} − (μ_m/σ_t) ∂(φ^{(0)}/2)/∂r`), so a closure
exact for an affine-in-`μ` flux is exactly the closure that gets the first-order
diffusion limit right. **⟹ BOTH of the brief's candidate phrasings are correct
statements of what the sources claim; "exact for a flux linear in μ" is the pointwise
form and "first-order diffusion-limit consistent" is its consequence.**

⚠ **What is NOT claimed:** BMC are explicit that τ is irrelevant at LEADING order.
Printed p. 159 (R-Z), immediately after Eq. (66): *"this result is true for any choice
of weighting factors in the approximation of the angular flux."* ⟹ the τ choice buys
you the FIRST-order diffusion limit only. Any τ ∈ [0,1] gets leading order.

## 2.3 — ⭐⭐ (c) THE CRUX. **`SUPPORTS` — BMC's cylinder analysis is graded on "ψ affine in the RADIAL (and axial) cosine" ONLY. The azimuthal cosine never enters.**

`[VERIFIED-ON-SCAN PDF pp. 11–13 = printed pp. 158–160]`

This is answered directly and unambiguously, in three printed places:

**(1) The direction vector they use has NO azimuthal component.** Printed p. 158,
**Eqs. (61)–(62)** — the `O(1)` relation that *defines* the first-order angular flux
in R-Z:

```
Ω⃗_{m,n} · ∇⃗ (φ^{(0)}/4π)  +  σ_t ψ^{(1)}_{m,n}  =  (1/4π) σ_t φ^{(1)}          (61)
Ω⃗_{m,n} = μ_{m,n} e_r  +  ξ_{m,n} e_z ;      ∇⃗ = (∂/∂r) e_r + (∂/∂z) e_z       (62)
```

⟹ solving (61), `ψ^{(1)}_{m,n} = (1/4π)φ^{(1)} − (1/4πσ_t)[μ_{m,n} ∂_r φ^{(0)} +
ξ_{m,n} ∂_z φ^{(0)}]`: **affine in `(μ, ξ)` = (RADIAL, AXIAL) cosines. The azimuthal
cosine `η` is absent from Eq. (62) entirely** (axisymmetry ⟹ `∇` has no `e_φ`
component ⟹ `J_φ ≡ 0`). **In 1-D cylindrical (`r` only) this collapses to affine in
the radial cosine ALONE** — exactly the brief's hypothesis.

**(2) The Galerkin ansatz the closure is graded against is affine in the radial cosine
alone.** Eq. (1), §2.2(i) above: `ψ_m = (1/4π)φ + (3/4π) J_r μ_m`, `J_r` = *radial*
current. One first-order mode, and it is radial.

**(3) The two first-order conditions they DO write are on the radial and axial cosine
edges — and the axial one is trivially satisfied.** Printed p. 159, **Eq. (69)** is a
PAIR of conditions: `Σ_n Σ_m (μ+ξ)(α_{m+½}μ_{m+½} − α_{m−½}μ_{m−½}) = 0` **and**
`Σ_n Σ_m (μ+ξ)(α_{m+½}ξ_{m+½} − α_{m−½}ξ_{m−½}) = 0`. **Eq. (70)** kills the second
one identically, and their reason is a one-liner:

> *"It is easy to show that the second relationship in Eq. (69) is zero because
> `ξ_{m−1/2,n} = ξ_{m+1/2,n} = ξ_{m,n}`."*

i.e. **the axial cosine is CONSTANT on a level, so it has no cell structure and needs
no closure.** ⟹ the only surviving condition, **Eq. (73)/(75)**, is a sum over the
**radial-cosine** cell edges alone. **There is NO azimuthal-cosine condition anywhere
in the paper.**

⟹ **VERDICT for the brief's crux question: `SUPPORTS`.** The literature grades a
cylindrical angular closure on its accuracy for **`ψ` affine in the radial cosine
`η`** (plus a trivially-satisfied axial mode in 2-D R-Z), **NOT** for `ψ` affine in
`(η, ξ)`.

**Does any source treat the azimuthal cosine `ξ` as an independent first-order mode?
`REFUTES` — no, and two sources state the symmetry that forbids it:**

- **Hébert 2009 Eq. (3.153)**, printed p. 91 `[VERIFIED-ON-SCAN]`:
  `φ(ρ,ξ,ω) = φ(ρ,ξ,−ω)` (together with `φ(ρ,ξ,ω) = φ(ρ,−ξ,ω)` in HIS letters, i.e.
  evenness in the axial cosine). `ω → −ω` flips `sin ω` and preserves `cos ω`, so in
  the brief's letters this is **exactly `ψ(η, ξ) = ψ(η, −ξ)`** — the reflection
  symmetry named in the brief. Its immediate consequence, which Hébert then uses:
  the source moments with `m < 0` vanish (3.154).
- **Carlson & Lathrop 1965, LA-3251, printed p. ~28** `[sidecar]`: *"In one-dimensional
  cylindrical geometry the flux is symmetric in `ξ` and `η`, so that only a quadrant,
  containing the entire range of `μ`, is needed."* In LA-3251's cylindrical letters
  `μ` = radial, `η` = **azimuthal**, `ξ` = axial ⟹ **the flux is EVEN in the azimuthal
  cosine, and the full range is needed only in the RADIAL cosine.** This is the
  published statement of the brief's premise, in the original S_N report.

⟹ `J_φ = ∫ ξ ψ d²Ω = 0` identically (odd integrand against an even flux). An
azimuthal first-order mode does not exist to be closed. **No local source treats `ξ`
as an independent first-order mode.**

## 2.4 — (d) THE CLAMP. **`REFUTES` a `τ ∈ [½,1]` limiter, three separate ways.**

**(i) BMC state the admissible range, and it is `[0,1]`.** Printed p. 157, immediately
under Eq. (53) `[VERIFIED-ON-SCAN PDF p. 10]`:

> *"where the weighting factors can take on any value between zero and one."*

followed by *"τ_{m,n} = 1 gives the step scheme and τ_{m,n} = ½ gives the diamond
scheme"* — so `½` is an *interior* landmark, not a floor. The identical sentence
appears for the sphere at Eq. (15).

**(ii) BMC's own worked S₂ example produces `τ₁ < ½`.** Printed p. 155,
**Eq. (47)** `[VERIFIED-ON-SCAN PDF p. 8 — see §2.5]`, obtained after Eq. (46) forces
`μ_{3/2} = 0`:

```
τ_1 = μ_1 + 1        and        τ_2 = μ_2                                      (47)
```

With the S₂ Gauss–Legendre ordinates `μ_1 = −1/√3`, `μ_2 = +1/√3`:
**`τ_1 = 1 − 1/√3 ≈ 0.42265` and `τ_2 = 1/√3 ≈ 0.57735`.**
⚠ **Provenance precision:** the paper prints the *formula* `τ_1 = μ_1 + 1`; the
*number* `0.42265` is one substitution away and is NOT printed. Cite it as
"BMC Eq. (47) with S₂ Gauss–Legendre", not as "BMC's value".
⟹ ORPHEUS's `max(0.5, min(1.0, τ))` floor would **clamp BMC's own published S₂
example**, re-introducing the β contamination τ exists to remove.

**(iii) BMC say what `τ ≡ ½` costs, and it is the quadrature.** Printed p. 155, right
after Eq. (47):

> *"A simple average relationship, where all τ values equal ½, would force the
> quadrature set to be the midpoint rule, μ₁ = −½, μ₂ = ½. However, the midpoint rule
> is not a good choice of quadrature sets in general."*

⟹ `τ ≡ ½` is not "the safe default"; it is a *constraint on the ordinates*. With any
other quadrature, `τ ≡ ½` breaks Eq. (43) and hence Eq. (41)/(75).

**No source in the corpus prescribes any limiter, clamp, floor or ceiling on τ.** The
only bound stated anywhere is `[0,1]`, and Lathrop 2000 gives it as a *consequence*,
not an imposition — printed p. 245, on his `δ` (with `τ = (1+δ)/2`): *"|δ| ≤ 1 as long
as `μ_m` is in `Δμ`"*, i.e. the bound is the statement that the ordinate lies inside
its own angular cell. The only negative-flux fixups in the corpus are **spatial**
(Hébert 3.387–3.389).

## 2.5 — ⭐⭐ THE STRUCTURAL FINDING the brief did not ask for but needs: **BMC's cylinder face carries TWO different cosines, and only one of them is the streaming coefficient**

This reconciles Claim 1(b) with Claim 2(a), which otherwise look contradictory ("the
face coefficient is `ξ`" vs "τ is barycentric in `η`"). Both are true because
**each azimuthal cell face carries two independent numbers**:

| at face `m±½` on level `n` | symbol | defined by | used for |
|---|---|---|---|
| **azimuthal** cosine (brief `ξ`) | BMC `α_{m±½,n}` (= `W_p η_{p,q±½}` in Hébert) | conservation recursion: BMC **(50)** / Hébert **(3.398)** / Carlson–Lathrop **(3-10)** | the **streaming/redistribution coefficient** — Claim 1(b) |
| **radial** cosine (brief `η`) | BMC `μ_{m±½,n}` | cumulative-weight recursion: BMC **(52)** | **defining τ** (Eq. 74) and evaluating β (Eq. 75) — Claim 2(a) |

⟹ **Claim 1(b) and Claim 2(a) refer to different quantities at the same face.** There
is no conflict. `[VERIFIED-ON-SCAN: BMC pp. 9–10, 13; Hébert p. 138]`

⚠ **And BMC (52) carries a HIDDEN QUADRATURE PREMISE that ORPHEUS's rule violates.**
Eq. (52) marches the radial cosine by the *weight*: `μ_{m+½,n} − μ_{m−½,n} = w̄_{m,n}`.
That is only the radial-cosine cell measure if **the within-level weight equals the
cell's radial-cosine width** — true by construction in BMC's level-symmetric family,
FALSE for an equispaced-ω rule. `[M]` on Hébert's own cylinder quadrature (see §2.6)
the arc cell's radial-cosine measure is `2 sinθ sin(ω_m) sin(Δω/2)` ∝ `sin ω_m`, while
the weight is CONSTANT — mismatched by up to ≈5× at M = 8. Feeding BMC (52) an
equispaced-ω rule therefore produces cells that **do not contain their own nodes**
(τ outside `[0,1]`, worsening with refinement). **Eq. (52) is not a law; it is the
statement that in BMC's quadrature family the weight equals the cell's radial-cosine
measure.** Eq. (74)/(43) itself is silent about edge provenance, so it survives a
different partition — Eq. (52) does not.

## 2.6 — ⭐ Hébert's OWN cylinder quadrature is equispaced-in-ω with EQUAL weights — the same family ORPHEUS uses — and he pairs it with recursion-defined edges + plain diamond

`[VERIFIED-ON-SCAN PDF p. 69 = printed p. 135]`, Eqs. (3.369)–(3.373):

- **(3.369)** `μ_{p,q} = √(1−ξ_p²) cos ω_{p,q}` — `ξ_p` = positive Gauss–Legendre base points.
- **(3.370)** `ω_{p,q} = (π/2)[1 − (N−2q−2p+3)/(N−2p+2)]` ⟹ `[M]` **cell-centred
  equispaced ω** on each level, `Δω = π/(N−2p+2)`, and **ω INCREASES with `q`**.
- **(3.371)** `W_{p,q} = π W_p /(N−2p+2)` ⟹ **all within-level weights EQUAL**, and
  **(3.372)** `Σ_p W_p = 1`, `Σ_p Σ_q W_{p,q} = π`.
- **(3.373)** the `N/2` zero-weight starting points sit at **`ω = π`**, with
  `μ̃_p ≡ μ_{p,N−2p+5/2} = −√(1−ξ_p²)` — the TOP half-index. Confirms §1.5(3).
- The alternative **product** quadrature (3.374)–(3.378) is the same structure.

`[M]` **On this quadrature Hébert's own edge recursion (3.398) reduces to
`η_{q+½} = η_{q−½} + Δω · √(1−ξ_p²) cos ω_q`** (because `W_{p,q}/W_p = Δω`), and its
output is **exactly `κ ×` the geometric arc value**, `κ = Δω/(2 sin(Δω/2))`, the SAME
`κ` at every interior face, closing exactly at `0` at the far end:

```
M = 8:  eta_recursion / (sinθ · sin ω_{q+1/2})  =  1.0064545428  at all 7 interior faces
        predicted κ = Δω/(2 sin(Δω/2))          =  1.0064545428      (agreement 1e-12)
        η at the ω = π end                       =  1.1e-16          (closes)
```

⟹ **`κ` is not an error in α.** It is the `O(Δω²)` gap between the *conservative*
face cosine (what α is, in all three sources) and the *geometric* arc face cosine.
Driving `κ → 1` would break conservation. **And Hébert's published angular closure on
this quadrature is the PLAIN DIAMOND, (3.406) — `τ ≡ ½`, no τ, no Morel–Montry.**
So "Hébert's scheme" and "BMC's scheme" are two complete, self-consistent, DIFFERENT
schemes; ORPHEUS mixes Hébert's quadrature with BMC's τ, which is legitimate *only*
via Eq. (43)/(74)'s edge-agnostic form, **not** via Eq. (52).

## 2.7 — ⚠⚠ THE SIGN TRAP in the closed-form geometric-edge τ: it is set by the MARCH ORIENTATION

`[M]` If the edges are taken at the geometric arc boundaries `ω_m ± Δω/2`, then Eq. (74)
gives a closed form whose sign depends entirely on which edge is labelled `m−½`:

```
index marches WITH ω      (Hébert's q, Eq. 3.370)   →  τ_m = ½ − ½ cot(ω_m) tan(Δω/4)
index marches WITH the RADIAL COSINE, i.e. AGAINST ω
       (BMC Eq. 52's μ_{1/2}=−√(1−ξ²) → +√(1−ξ²);  = ORPHEUS)  →  τ_m = ½ + ½ cot(ω_m) tan(Δω/4)
```

`[M]` both verified to `8.9e-16` at M = 8, `sinθ = 0.7`. **The two produce the SAME SET
of τ values in REVERSED ORDER** (`0.2524 … 0.7476` vs `0.7476 … 0.2524`), and both
satisfy `τ ∈ [0,1]` at every M. ⟹ **ORPHEUS's `+` sign is the correct one for its
η-ascending storage** (BMC's ordering). A port of Hébert's `q` indexing would need the
`−` sign; getting this backwards silently reverses τ across the level, which is
exactly the kind of error a symmetric fixture cannot see.

## 2.8 — Adjudication of the PRIOR deliverable's line ~943 claim

The sentence under review, from `scratch/q64_tau_edge_convention_literature.md` §F4:

> *"`κ` is … the ratio between **the recursion-defined edge — which is what α is, by
> definition, in both sources — and the geometric arc-half-angle edge, which no source
> uses.**"*

**It is TWO claims, and they grade differently.**

**(i) "the recursion-defined edge is what α IS, by definition, in both sources" —
`SUPPORTS`, and now in THREE sources, including the original.**
Hébert **(3.398)** defines `η_{p,q±½}` by a recursion from constant-flux preservation;
BMC **(50)** defines `α_{m±½,n}` by a recursion from the same requirement;
Carlson & Lathrop **LA-3251 (3-10)** does it first, from "everywhere uniform flow".
`[M]` **no equation in Hébert §3.9.3 ever computes an azimuthal half-angle
`ω_{p,q±½}`** — a grep of the whole chapter finds that symbol only inside the flux's
*argument* in Eq. (3.393), i.e. as a label for *which direction* the half-angle flux
belongs to. `[VERIFIED-ON-SCAN p. 72]`

**(ii) "the geometric arc-half-angle edge, which no source uses" — this is the PRIOR
AGENT'S OWN INFERENCE, not a sourced statement, and its quantifier is unsafe.**
- No source *says* it. It is an absence-of-evidence claim.
- Its denominator was **2** (Hébert §3.9.3 + BMC's S₂ example). I have now widened it to
  **3 of 3** local sources that define a face cosine — LA-3251 also uses a recursion —
  and I found **no counterexample** anywhere in the local corpus (grep for half-angle /
  bisection / `ω_{m±½}` returns only Hébert's flux labels and Stacey's unrelated
  *half-angle Legendre* `D-P_L` polynomials).
- ⛔ **But two PRIMARY sources that would settle it are NOT LOCAL**, and both are named
  by the local text as the actual authorities:
  - **Alcouffe & O'Dell** — Hébert's ref. **[36]**, credited by name at printed p. 138
    (*"We now present the approach proposed by Alcouffe and O'Dell to compute the
    `η_{p,q±1/2}` values"*). **This is the primary source for the cylinder edge
    construction and it has never been read.** Hébert's Ch.-3 reference list is not in
    the local PDF (chapter extract only), so even the full citation is unresolved.
  - **Morel & Montry (1984), TTSP 13(5):615** — the origin of τ itself, credited
    throughout BMC. Also not local.
- ⚠ **A sharper statement is available and is what I would put in a theory page**, since
  it is a claim about *mechanism* rather than a survey of all literature:

  > **No source in the corpus evaluates a cosine at an azimuthal half-ANGLE. All three
  > that define a face cosine define it by a conservation or cumulative-weight
  > recursion. Whether that recursion's output coincides with the geometric arc value
  > is a property of the QUADRATURE, not of the scheme** — for BMC's level-symmetric
  > family the radial-cosine march (52) lands on it exactly; for an equispaced-ω rule
  > Hébert's azimuthal-cosine march (3.398) lands on `κ ×` it.

  ⟹ the prior deliverable's "which no source uses" is **not wrong on the evidence but
  is stated at the wrong quantifier**, and it obscures that BMC's Eq. (52) edges *are*
  the geometric radial-cosine edges whenever the quadrature makes weight = cell
  measure. Restate it as above.

## 2.9 — CLAIM 2 SUMMARY

| sub-claim | verdict | citation |
|---|---|---|
| (a) a cylinder τ formula exists | **`SUPPORTS`** — BMC **Eq. (74)**, printed p. 160 | edges **(52)**, closure **(53)** |
| (a) barycentric in the RADIAL cosine `η`, not in `ω` | **`SUPPORTS`** ORPHEUS; **`REFUTES`** an ω-barycentric reading | BMC `μ = sinθ cos ω` (p. 156) |
| (b) "exact for a flux linear in the radial cosine" | **`SUPPORTS`** | BMC sentence under **(43)**, with **Eq. (1)** + **(15)** |
| (b) "first-order diffusion-limit / Galerkin consistency" | **`SUPPORTS`** — the same property, other side | BMC after **(74)**; **(75)**; uniqueness claim p. 161 |
| (b) leading order needs no τ | **`SUPPORTS`** | BMC p. 159 after **(66)** |
| (c) cylinder diffusion-limit flux affine in the RADIAL cosine only (axial mode trivial, azimuthal absent) | **`SUPPORTS`** — the crux, answered 3 ways | BMC **(61)/(62)**, **(69)/(70)**, **(1)** |
| (c) any source treating azimuthal `ξ` as an independent first-order mode | **`REFUTES` — none** | Hébert **(3.153)**; LA-3251 printed p. 31 |
| (d) a limiter/clamp such as `τ ∈ [½,1]` | **`REFUTES` — no source has one; range is `[0,1]`** | BMC under **(53)** and **(15)**; Lathrop p. 245 |
| (d) BMC's own S₂ value | **`τ₁ = 1 − 1/√3 ≈ 0.42265 < ½`** — below any `½` floor | BMC **Eq. (47)** + **Fig. 1** ordinates |

---

# NOTATION CLASHES — the consolidated flag list

| # | trap | consequence if missed |
|---|---|---|
| N1 | **Hébert / BMC / LA-3251 (cyl): `μ` = RADIAL, `η` = AZIMUTHAL, `ξ` = AXIAL.** Brief/ORPHEUS: `η` = radial, `ξ` = azimuthal, `μ_z` = axial. | Every coefficient in Claim 1 reads inverted. |
| N2 | **B&G 1970 and Stacey 2007: `μ` = AXIAL cosine**, azimuth is `χ` (B&G) / `φ` (Stacey). B&G reserve `ω` for the *spherical* azimuth. | `√(1−μ²)cosχ` mistaken for a polar factor. |
| N3 | ⛔ **BMC printed p. 156 says `η = sinθ cos ω` — a PUBLISHED TYPO** (must be `sin ω`). | The azimuthal-face coefficient read as the radial cosine — i.e. Claim 1 read backwards, from the paper ORPHEUS cites. |
| N4 | ⛔ **Stacey Eq. (9.150)'s `sinθ/4` — a PUBLISHED TYPO** (must be `sinφ/r`); dimensionally impossible as printed. | A "second source" that appears to contradict Claim 1. |
| N5 | ⛔ **BMC (50) and (52) are printed self-referentially** (`α_{m+½}=α_{m+½}−…`, `μ_{m+½}=μ_{m+½}+…`); RHS must be `m−½`. Sphere twins (11)/(12) are printed correctly. | A vacuous recursion transcribed into a solver. |
| N6 | **`β` names THREE different objects**: BMC's scalar `J^{(2)}` contamination coefficient (41/75); Lathrop-2000's *sequence* = α's pointwise defect (25); **LA-3251's `β` = the y-curvature redistribution coefficient (3-18)**. | "β = 0" asserted about the wrong object. |
| N7 | **α SIGN is not universal**: BMC (50) and LA-3251 (3-10) use `−w μ`; Hébert (3.399) uses `+W μ`. Tracks each author's sign in front of the redistribution term. | α built with the wrong sign against ORPHEUS's own streaming sign. |
| N8 | **Cylinder redistribution sign is `−`, sphere's is `+`** — physical (ω decreases along the path; μ increases), not conventional. | A sign ported from the sphere. |
| N9 | **The geometric-edge τ closed form flips sign with the march orientation** (§2.7). | τ silently reversed across the level. |
| N10 | **LA-3251's letters PERMUTE between geometries** (footnote, printed p. 17): the cylindrical assignment is not the spherical one. | Cross-geometry mis-mapping inside one report. |

---

# PROVENANCE / TOOLING NOTES

- **Zotero DOWN this session.** `nc -z 127.0.0.1 23119` → **refused**. Per
  `.claude/rules/delegation.md` + L-007 this is a broken server, not an empty library.
  **No user annotations were checked.** If the τ convention is going into a theory page,
  re-query Zotero when it is up — user highlights on BMC Eq. (74)/(52) would be
  high-signal on exactly the §2.5 question.
- **OSTI unreachable** this session (`ConnectionResetError`); CrossRef and OpenAlex
  worked. Status of the two non-local primaries:
  1. ✅ **Morel & Montry 1984 — RESOLVED and now VERIFIED against CrossRef.**
     J. E. **Morel** and G. R. **Montry**, *"Analysis and elimination of the
     discrete-ordinates flux dip"*, **Transport Theory and Statistical Physics 13(5),
     615–633 (1984)**, **DOI `10.1080/00411458408211661`**, cited-by 19. ⟹ the
     volume/issue/page carried in prior notes was CORRECT; it is no longer an
     unverified inheritance. **This is the origin paper for τ and for the β = 0
     condition, and it is the one document that would settle the clamp question at the
     primary source. NOT in `scratch/literature/` — worth acquiring.**
  2. ⛔ **Alcouffe & O'Dell — UNRESOLVED. Do NOT cite.** Hébert's ref. **[36]**, named
     at his printed p. 138 as the authority for the cylinder `η_{p,q±½}` construction.
     Four searches (OSTI ×2 unreachable, OpenAlex ×2, CrossRef ×2, INIS ×1) returned no
     match. The local PDF is a chapter extract with no reference list, so the citation
     string is not recoverable here either. **Best lead, NOT confirmed:** CrossRef does
     confirm a *"Handbook of Nuclear Reactors Calculations"* (Ronen, ed.) reviewed in
     TTSP 1988 — a plausible venue for an Alcouffe/O'Dell chapter, but **there is no
     evidence they authored one**, so treat any such attribution as a hypothesis.
     Resolving it needs Hébert's full book bibliography (p. ~189 of Ch. 3, not in the
     local extract). ⚠ Until resolved, **the primary source for ORPHEUS's cylinder
     angular-cell-edge construction has never been read by anyone on this project.**
- ⚠ **Provenance correction for ORPHEUS's code comment.** The cylinder τ should cite
  **BMC Eq. (74)** (with edges (52), closure (53)); "BMC Eq. 43 = Lathrop Eq. 23" is the
  **spherical** pair. Same content, wrong geometry's equation number.
