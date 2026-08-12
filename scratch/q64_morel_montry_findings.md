# Morel & Montry (1984) — the PRIMARY source for the angular weighted diamond

**Citation** (verified against the article's own title page, PDF p. 2):
J. E. Morel and G. R. Montry, "Analysis and Elimination of the Discrete-Ordinates
Flux Dip," *Transport Theory and Statistical Physics* **13**(5), 615–633 (1984).
DOI `10.1080/00411458408211661`. Sandia National Laboratories. A revision of a
paper given at the ANS topical "Advances in Reactor Computations," Salt Lake
City, March 1983 (footnote, PDF p. 2 / printed 615).

- PDF: `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/literature/Morel-Montry(1984)Analysis and elimination of the discrete-ordinates flux dip.pdf` (20 pp)
- Sidecar: `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/literature_ocr/Morel-Montry(1984)Analysis and elimination of the discrete-ordinates flux dip.md`
- **Page map: printed p. = PDF p. + 613** (`[M]` verified at both ends: PDF p. 2
  header = printed 615, PDF p. 20 header = printed 633). Cited below as
  `PDF p. N / pr. N+613`.

### Denominator for every "NOT ADDRESSED" verdict

The **entire** paper was read. 19 printed pages, structure: Abstract ·
Introduction · "Analysis of the Flux Dip" · "Elimination of the Flux Dip" ·
Summary · 3 References · Appendix (cylinder). PDF pp. 9, 14, 15 are figure-only;
Figs. 1–8 + A1, A2 are bitmaps with captions carried in the sidecar. **There are
no tables** (`.mocr.json` carries zero `tables` entries — L-008 checked). The
complete equation inventory is (1)–(19) with sub-numbers (1a)–(1d), (4a)–(4b),
(6a), (7a), (16a)–(16b), (17a)–(17b), plus (A1)–(A4). **Nothing else exists in
which a clamp, an endpoint rule, or a positivity condition could hide.**

### Pages verified visually against the scan (SSOT)

PDF **4** (pr. 617, the accuracy claim + Eqs. 1/1a) · **11** (pr. 624, Eqs.
12/13/14 + the seed threshold) · **12** (pr. 625, Eqs. 15/16a/16b/17a/17b —
the sphere τ) · **17** (pr. 630, Summary end + References + Appendix opening) ·
**20** (pr. 633, Eqs. A1–A4 — the cylinder τ). The sidecar was faithful on all
five except where noted under "Errata" at the end.

---

## Q1 — Does M&M prescribe ANY limiter, clamp, bound, or special case on τ?

**Source.** Sphere: PDF p. 12 / pr. 625, Eqs. (16a)–(16b). Cylinder: PDF p. 20 /
pr. 633, Eqs. (A1)–(A2). Both **visually verified**.

τ is introduced and left bare. Eq. (16b) defines it as a ratio of two cosine
differences and the text moves straight on to proving β = 0; Eq. (A2) repeats
the identical ratio for the cylinder and the paper then **ends** (A4 is the last
equation on the last page). There is no interval statement, no sign condition,
no "if the cell is degenerate" branch, no fallback to plain diamond or step, no
mention of any τ value being inadmissible, and no numerical-stability remark
anywhere in the 19 printed pages.

The one *structural* fact worth recording in place of a clamp: because the
edge march (15) is monotone in μ and the quadrature node lies inside its own
cell, τ from (16b) is **automatically in (0,1)** — a clamp to [0,1] would be
vacuous, and no *tighter* interval is ever proposed. `[M]` for the cylinder the
same holds: cos is monotone on the marched azimuthal range, so (A2) ∈ (0,1) by
construction.

`VERDICT: NOT ADDRESSED.` No limiter, clamp, bound, interval, sign condition, or
degenerate-cell special case exists anywhere in Morel & Montry (1984). Combined
with the prior reads this makes the clamp denominator **5 of 5 primaries, none
of which clamps τ** (M&M 1984, Reed–Lathrop 1970, Bailey–Morel–Chang 2010,
Hébert 2009, Lathrop 2000). The `[½,1]` interval remains Grant's, on the
*spatial* weight.

---

## Q2 — The endpoint / grazing ordinates; τ in a narrowing cell

**Source.** PDF pp. 3–4 / pr. 616–617 (Miller–Alcouffe); PDF pp. 8, 10–11 /
pr. 621, 623–624 (Eqs. 10–14, the effective starting cosine); PDF p. 12 /
pr. 625 (Eq. 17b, δ_m).

M&M discuss "the first direction" at length, but they mean the **starting
direction** μ_{1/2} = −1 — the *edge* that carries no weight and is solved with
the slab equation — not the first weighted **ordinate**. Their entire endpoint
apparatus is about that: the Miller–Alcouffe procedure, the "effective starting
cosine" μ_s of Eq. (10), Fig. 4's plot of μ_s(r), the inconsistency Eqs.
(11)–(13), and the seed generalisation Eq. (14). Nothing anywhere treats the
first or last *weighted* ordinate m = 1 / m = N as a special case of the τ
formula.

Nor is there any narrow-cell analysis. Eq. (17b) defines
δ_m = μ_m − (μ_{m+1/2}+μ_{m−1/2})/2 — precisely the node's offset from its own
cell centre, the quantity whose smallness R&L Eqs. (15)/(16) require — and M&M
use it **only** as an algebraic intermediate inside the parity proof of β = 0.
The single sentence in the paper that comments on its size is a flat concession
(PDF p. 12 / pr. 625): with standard sets the quadrature points are *never* at
the cell centres. They do not bound δ_m, do not relate it to the weight, and do
not ask what happens when the cell measure shrinks faster than the weight.

⭐ `[M]` **And the effect is already visible in their own sphere.** Reproducing
Eq. (16b) on the Gauss sets they use throughout (`sum(W)=1`, edges from Eq. 15):

| N | τ₁ | τ_N | τ at mid-level | (τ₁−½)/(2W₁) |
|---|---|---|---|---|
| 2 | 0.422650 | 0.577350 | — | −0.077 |
| 8 | 0.392282 | 0.607718 | 0.505770 | −1.064 |
| 32 | 0.389840 | 0.610160 | 0.500390 | −15.695 |
| 64 | 0.389708 | 0.610292 | 0.500099 | −61.848 |

The interior converges to ½ but **the end ordinates do not** — τ₁ saturates at
≈ 0.3897 and the R&L ratio |τ−½|/w **diverges as O(N)** at the ends, on the
sphere, with the plain Gauss quadrature M&M themselves use. So the
endpoint conflict is **not** a cylinder artefact of the σ_y fold; it is generic
to barycentric-in-μ τ wherever the end cells are weight-starved, and it is in
the primary source's own default configuration.

`VERDICT: NOT ADDRESSED (endpoint ordinates), and the phenomenon is present but
unremarked in their own sphere results.` They discuss the starting *direction*
exhaustively and the first weighted *ordinate* not at all.

---

## Q3 — The CYLINDER: what does M&M itself do?

**Source.** Appendix, PDF pp. 17–20 / pr. 630–633. Eqs. (A1)–(A4) **visually
verified on PDF p. 20**.

**Yes, they treat it — and in more detail than "analogous to".** The appendix is
the extension of both the Miller–Alcouffe procedure and the weighted diamond to
1-D cylindrical geometry, and its central paragraph (PDF p. 19 / pr. 632) states
the structural difference explicitly: in cylindrical geometry the cell-edge
**cosines** are *not* defined by partial sums of the weights as they are in the
sphere; it is the cell-edge **azimuthal coordinates** that are so defined, and
the edge cosines follow from the cell-centre polar and cell-edge azimuthal
coordinates. Their reason is stated in one clause: the angular derivative is
taken with respect to the azimuth in the cylinder and with respect to the cosine
in the sphere.

The scheme, per ξ-level ℓ, with m = 1 the weighted direction nearest the
starting direction:

- **(A1)** `ψ_{ℓ,m} = (1−τ)ψ_{ℓ,m−1/2} + τ ψ_{ℓ,m+1/2}`
- **(A2)** `τ = (μ_{ℓ,m} − μ_{ℓ,m−1/2}) / (μ_{ℓ,m+1/2} − μ_{ℓ,m−1/2})`
- **(A3)** `μ_{ℓ,m+1/2} = sin(θ_ℓ) cos(φ_{ℓ,m+1/2})`
- **(A4)** `φ_{ℓ,m+1/2} = φ_{ℓ,m−1/2} + π C_ℓ W_{ℓ,m}`, `φ_{ℓ,1/2} = −π`,
  where `C_ℓ` normalises the level's weights to unity.

**Answers to the three sub-questions:**

1. **τ is posed in the RADIAL COSINE μ** (A2), not in the azimuth. It is the
   same barycentric ratio as the sphere's (16b), letter for letter.
2. **The PARTITION is in the AZIMUTH** (A4): cell-edge azimuths march from −π by
   weight-proportional increments; the edge *cosines* are the cos-image (A3) of
   those azimuthal edges. So the cylinder has a two-stage construction —
   partition in φ, barycentre in μ = sin θ cos φ.
3. **Quadrature assumed:** whatever set the user brings, *provided* the
   azimuthal cell width is proportional to the weight (A4). **Their derivation
   does NOT depend on weight = cell μ-measure** in the cylinder; it depends on
   **weight = cell AZIMUTHAL measure** (per-level, normalised by C_ℓ). In the
   sphere the dependence is the other one — Eq. (15), weight = ½ × cell
   μ-measure, credited to Bell & Glasstone (ref. 3).

⭐⭐ **This directly ratifies the ORPHEUS 6.4 partition.** For an equal-weight
level (`W_{ℓ,m} = 1/M`), (A4) collapses to **equispaced azimuthal edges of width
Δφ = π/M with the node at the arc midpoint** — i.e. exactly the "midpoint in ω"
partition — and (A2)+(A3) then *determine* τ. `[M]` the closed form implied by
(A2)–(A4) on such a level is

```
τ_m = ½ − ½ cot(φ_m) tan(Δφ/4)          (φ marched −π → 0, node at the arc midpoint)
```

reproduced to machine precision against a direct evaluation of (A2)–(A4) at
M = 2, 4, 8, 16, 32, 64. It is **independent of the level's sin θ_ℓ** (the polar
sine cancels out of the ratio) — the same τ on every ξ-level. This is ORPHEUS's
`τ = ½ ± ½cot(ω)tan(Δω/4)`; the sign is purely the march orientation (M&M march
from the grazing direction at φ = −π toward φ = 0, so *their* τ₁ ≈ ¼ and
τ_M ≈ ¾).

⭐ **And it produces the ¼ endpoint, refinement-independently.** `[M]`:

| M | τ₁ | τ₂ | τ mid | τ_M | (τ₁−½)/w |
|---|---|---|---|---|---|
| 4 | 0.2599 | 0.4588 | 0.5412 | 0.7401 | −1.00 |
| 8 | 0.2524 | 0.4263 | 0.5098 | 0.7476 | −1.98 |
| 32 | 0.25015 | 0.41725 | 0.50060 | 0.74985 | −7.995 |
| 128 | 0.250009 | — | — | 0.749991 | −31.999 |

τ₁ → ¼ and τ_M → ¾ exactly, and (τ₁−½)/w = −M/4 — the O(M) divergence the brief
measured. `[M]` the end cell's μ-width relative to a mid-arc cell shrinks like
1/M² (0.1989 at M=8 → 0.01227 at M=128), confirming quadratically-narrow-in-μ
against linear weight. **M&M's own formula generates exactly ORPHEUS's observed
behaviour, and they neither flag it, bound it, nor mention it.**

⚠ Two further cylinder facts worth recording:

- **There is NO cylinder β analysis.** The appendix gives the recipe and nothing
  else — no diffusion-limit expansion, no Eq.-(6a) analogue, no proof that the
  cylinder scheme zeroes anything. The only cylinder evidence in the paper is
  empirical: they tested it and the dip was similarly eliminated (PDF p. 17 /
  pr. 630). Every "β = 0" claim in this paper is **spherical**.
- **The cylinder Miller–Alcouffe rule** (PDF pp. 18–19 / pr. 631–632): the flux
  at the cylinder centre must be azimuth-independent but may depend on polar
  angle, so all origin fluxes on a ξ-level are set to that level's starting-flux
  value, and each starting direction uses the same slab equation as in the
  sphere.

`VERDICT: TREATED, IN ITS OWN APPENDIX. τ is barycentric in the RADIAL COSINE
(A2); the cell PARTITION is in the AZIMUTH by cumulative weight (A4) and mapped
through cos (A3). The derivation depends on weight = cell azimuthal measure, NOT
weight = cell μ-measure. On an equal-weight level this IS the equispaced-ω
midpoint partition, and it yields τ₁ → ¼ / τ_M → ¾.`

---

## Q4 — The derivation of β = 0: what it constrains, and whether τ is unique

**Source.** PDF pp. 4–6 / pr. 617–619 (Eqs. 1–7a); PDF pp. 12–13 / pr. 625–626
(Eqs. 15–19).

**It is not a Galerkin derivation.** The trial space is a single-term
diffusion-limit ansatz, Eq. (2): `ψ_m = φ + 3J μ_m` — linear in μ, with φ and J
the scalar flux and current. It is admissible only if the quadrature meets the
three standard requirements they list at PDF p. 5 / pr. 618: cosines symmetric
about μ = 0, `Σ W_m = 1`, and the **diffusion condition** Eq. (3)
`Σ_m μ_m² W_m = 1/3`. Substituting (2) into the angle-discretised S_N equation
(1) gives (4); the **zeroth** quadrature moment gives the exact continuity
equation (5), and the **first** moment gives Fick's law with one spurious term,
Eq. (6), whose coefficient is

- **(6a)** `β = 3 Σ_{m=1}^{N} μ_m ( α_{m+1/2} μ_{m+1/2} − α_{m−1/2} μ_{m−1/2} )`

with `μ_{m±1/2}` the **scheme-implied** edge cosines — Eq. (4a) for diamond
(`μ_{m+1/2} = 2μ_m − μ_{m−1/2}`, seed `μ_{1/2} = −1`), Eq. (4b) for step
(`μ_{m+1/2} = μ_m`). This confirms the L-011 reading: (6a)'s half-index symbols
are the *scheme's* edges, not the quadrature's. β enters only through the
diffusion coefficient, Eq. (7a): `D = 1 / (3(σ_t + 2β/r))`.

**Which condition actually picks τ — and it is NOT β = 0.** The decisive
paragraph is PDF p. 12 / pr. 625: with diamond differencing the scheme-implied
edges coincide with the standard weight-defined edges (15) **only if every
quadrature point sits at its cell centre**; with step they never coincide; with
standard sets the points are never at the centres. M&M then choose the weighted
diamond precisely so that **coincidence is achieved for any quadrature set** —
that is, they impose the scheme-implied edges ≡ the standard edges (15), which
is **N equations, one per ordinate**, and inverting them ordinate-by-ordinate
*is* Eq. (16b). β = 0 is then **derived as a consequence**, Eqs. (17a)–(19), by a
pure parity argument: δ_m is odd, (α_{m+1/2}+α_{m−1/2}) is even, so every term in
(19) sums to zero over a symmetric level.

⭐⭐ **So BMC's framing is an inversion, and it is quantitatively wrong.** β is
**one scalar** for the whole quadrature; it cannot determine N unknowns.
`[M]` demonstrated on Gauss-S8 by root-finding τ₁ with the other seven τ's
randomised: three distinct τ-vectors with `β = O(1e-16)` and
`‖τ − τ_M&M‖_∞ = 0.238, 0.308, 0.242` respectively. β = 0 defines an
(N−1)-dimensional solution set; the barycentric τ is **one point on it**, picked
by the *stronger, per-ordinate* edge-coincidence requirement. (`[M]` β is a
genuine non-degenerate diagnostic in the sphere — random fold-symmetric τ gives
β ∈ [−0.30, +0.07], not round-off — so the sphere case is unlike the cylinder
equispaced-ω level where β is blind.)

`[M]` **The parity proof reproduces exactly.** Eq. (16b)'s τ gives β = 0.000e+00
at N = 2, 4, 8 and ≤1.6e-16 at N = 16, 32. Their other two claims reproduce too:
β_step > 0 and decaying with order (0.577 at S2 → 0.0427 at S32 — "β converges to
zero in the limit as the quadrature order is increased", PDF p. 7 / pr. 620), and
β_diamond = +0.1547 at Gauss-S2 (their Fig. 3 statement that "the β coefficient
is positive", PDF p. 8 / pr. 621), turning negative from S4 up (−2.68e-3 at S4,
−3.43e-5 at S8).

`VERDICT: β = 0 is a CONSEQUENCE, not the defining condition. The defining
condition is per-ordinate coincidence of the scheme-implied edges with the
weight-defined edges (15), which uniquely fixes each τ_m. β = 0 alone is one
scalar equation with an (N−1)-dimensional solution family — CONTRADICTS the
established list's paraphrase of BMC line 657.`

---

## Q5 — The starting direction and the α march

**Source.** PDF p. 4 / pr. 617, Eq. (1a) **visually verified**; PDF pp. 3–4 /
pr. 616–617 (Miller–Alcouffe); PDF pp. 8, 10–11 / pr. 621, 623–624 (Eqs. 10–14)
**p. 11 visually verified**; PDF p. 16 / pr. 629 (coarse-mesh truncation).

**The α march.** Eq. (1a), verified on the scan:
`α_{m+1/2} = α_{m−1/2} − μ_m W_m`, with **both** `α_{1/2} = 0` **and**
`α_{N+1/2} = 0` printed on the same line. Sign and form match ORPHEUS exactly.
They state both terminal conditions and say **nothing** about which end to march
from, nothing about accuracy or conditioning of the recurrence, and nothing about
error accumulation along it. (Note their weight normalisation is `Σ W_m = 1`, so
their `W_m` is half the `Σw = 2` convention's weight; the α values scale
accordingly.)

**The starting direction is a separate, non-α object.** The Miller–Alcouffe
procedure (PDF pp. 3–4 / pr. 616–617) is (i) use the **slab-geometry equation**
for the starting-direction flux, and (ii) set **all** weighted flux values at the
origin to that starting-flux value. M&M adopt it wholesale — every result in the
paper, including the dip-free ones, is "weighted diamond **in conjunction with**
the Miller–Alcouffe procedure", and the Summary (PDF p. 16 / pr. 629) attributes
the dip to a *three-way* interaction: the angular differencing scheme, the
calculation of the starting flux, and the boundary treatment at the origin.

⭐⭐ **The high-value part for a seed change: they quantify how a seed error
propagates.** Eq. (14) generalises the edge-march seed from −1 to the *effective*
starting cosine μ_s defined by Eq. (10) (`ψ̃_s = φ̃ + 3J̃ μ_s`, measured from an
actual S₂ run). Then, PDF p. 11 / pr. 624 **verified on the scan**: for a
Gauss-S₂ set with angular diamond differencing, **β becomes negative for
μ_s < −2/√3**.

`[M]` Reproduced exactly, and the dependence is *linear and un-damped*:

```
β(μ_s) = μ_s + 2/√3          (Gauss-S₂, angular diamond, Eq. 6a with Eq. 14)
```

verified at μ_s ∈ {−1.4, −1.2, −1.0, −0.8}; β(−2/√3) = −2.2e-16, β(−1) = +0.1547.
**The seed enters β with unit gain and no attenuation** — a seed perturbation
δμ_s moves β by δμ_s, and β is what sets the diffusion coefficient (7a) and hence
the dip. There is no averaging or damping along the march to soften it. Their
diagnostic instrument for this is the plot of μ_s(r) (Figs. 4, 6, 8); the
signature of a healthy calculation is μ_s transitioning to a constant −1 in the
interior (Fig. 6), and the S₂-diamond failure is μ_s dropping to ≈ −1.35 at the
origin (Fig. 4, PDF p. 10 / pr. 623).

**Two more seed-adjacent findings:**

- **The sign of the seed error decides whether you get a dip at all.** A seed
  that *underestimates* the starting flux (μ_s < −2/√3) drives β negative and
  produces the dip; one that *overestimates* it drives β positive, which
  destroys the exact diffusion-limit correspondence but **cannot** produce a dip
  (PDF pp. 7, 16 / pr. 620, 629). Their closing rule of thumb (PDF p. 17 /
  pr. 630) is that the dip is eliminated with any spatial scheme as long as the
  starting flux is not seriously underestimated relative to the weighted fluxes.
- **The starting direction has a DIFFERENT spatial truncation order** (PDF
  p. 16 / pr. 629, offered as a speculation): the spatial truncation error goes
  like Δr²/r² for the weighted fluxes but like Δr² for the starting flux — so on
  a coarse mesh the two are inconsistent, the imbalance grows toward small r, and
  β(r) becomes effectively non-zero (positive) even with the weighted diamond.
  They note the "which one is wrong" framing is a matter of convention, but that
  the starting flux is the *more* accurate of the two.

`VERDICT: BOTH terminal conditions are stated (α_{1/2} = 0 AND α_{N+1/2} = 0);
NO statement on march direction, conditioning, or recurrence error growth. But
the SEED of the angular-cell-edge march is analysed in depth: Eq. (14) admits an
arbitrary seed μ_s, and [M] β = μ_s + 2/√3 for Gauss-S₂ diamond — unit gain, no
damping, sign-determining for the dip. Their diagnostic is the effective starting
cosine, Eq. (10).`

---

## Q6 — Negative half-angle fluxes / positivity

**Source.** PDF pp. 6–7 / pr. 619–620 (the dip mechanism); whole paper searched
for a positivity statement.

The paper's negativity analysis is entirely about the **diffusion coefficient**,
never about the angular flux. The mechanism (PDF p. 7 / pr. 620) is precise and
worth having: if β < 0 and the point `r_D = −2β/σ_t` lies inside the sphere,
then D from (7a) is negative for r < r_D, **unbounded at r = r_D**, and positive
beyond. Since the current must be positive everywhere in a purely scattering
sphere (every entering particle must leave) and J = −D dφ/dr, the flux slope must
be positive inside r_D, zero at r_D, and negative outside — which *is* the dip.
`r_D` predicts the dip location, confirmed against the S₂ calculation to within
the spatial mesh resolution (PDF p. 11 / pr. 624).

So the "dip" is a **non-physical local minimum of the scalar flux**, not a
negative flux. There is no statement anywhere about ψ_{m±1/2} ≥ 0, no analysis of
when the angular recurrence produces negatives, no positivity-preserving
constraint on τ, and no fix-up scheme.

`VERDICT: NOT ADDRESSED. No positivity condition on the half-angle angular flux
exists in this paper. The only negativity analysed is of the S_N diffusion
coefficient D (Eq. 7a) via β < 0, which is the dip mechanism, not a flux-sign
condition. M&M provide NO authority for a half-angle positivity gate.`

---

## Q7 — Anything that CONTRADICTS the established list

**Confirmed by the primary (no contradiction):**

1. **τ IS barycentric-in-μ in the primary source** — Eq. (16b) for the sphere
   and (A2) for the cylinder, both visually verified. This is *not* a BMC
   re-derivation. BMC Eq. 42/43 reproduces M&M (16b) letter for letter; BMC
   Eq. 74 reproduces M&M (A2). Its SHAPE: numerator = (node cosine − lower edge
   cosine), denominator = (upper edge cosine − lower edge cosine), **in the
   direction cosine μ in both geometries**.
2. **The α recursion** `α_{m+1/2} = α_{m−1/2} − μ_m W_m` with `α_{1/2} = 0` is
   the primary's Eq. (1a), sign-for-sign as ORPHEUS ships it.
3. **The closure orientation** matches: M&M (16a) is
   `ψ_m = (1−τ)ψ_{m−1/2} + τψ_{m+1/2}`, i.e. τ weights the **upper/outgoing**
   edge, identical to ORPHEUS's `ψ_m = τ_m ψ_{m+1/2} + (1−τ_m)ψ_{m−1/2}`.
4. **The R&L midpoint criterion is acknowledged and knowingly violated.** M&M's
   δ_m (17b) *is* R&L's node-offset, and M&M state plainly that standard sets
   never put the node at the cell centre (PDF p. 12 / pr. 625).

**Contradictions / corrections to the established list:**

5. ⛔ **"Forcing β to zero determines the Morel and Montry weights" (BMC line
   657) is FALSE as stated.** See Q4: β is a single scalar, `[M]` shown to admit
   an (N−1)-dimensional family of τ-vectors. M&M derive τ from **edge
   coincidence** (N conditions) and prove β = 0 as a corollary. The established
   list should read "β = 0 is *implied by* the M&M weights", not "determines"
   them. BMC's β (their Eq. 41) is a first-order diffusion-limit *diagnostic*
   that the M&M τ zeroes — a necessary, not sufficient, condition.
6. ⚠ **Weight-normalisation mismatch.** M&M use `Σ W_m = 1` over μ ∈ [−1,1]
   (Eq. 1d) and the diffusion condition `Σ μ² W = 1/3` (Eq. 3), so their
   Eq. (15) reads `μ_{m+1/2} = μ_{m−1/2} + 2W_m` — **weight = HALF the cell
   μ-measure**. Any comparison against a `Σw = 2` convention must carry the
   factor 2, including in Eq. (1)'s `2/r` redistribution prefactor and in the α
   magnitudes.
7. ⚠ **The cylinder does NOT inherit "weight = cell μ-measure".** In the
   cylinder the weight is the cell **azimuthal** measure (A4); the μ-edges are
   the cos-image (A3). M&M say this in as many words (PDF p. 19 / pr. 632). Any
   statement of the form "M&M require weight = cell μ-measure" is a sphere-only
   statement.
8. ⭐ **M&M state their own order, and it is FIRST — deliberately.** PDF p. 4 /
   pr. 617, **visually verified**: unlike the Reed–Lathrop scheme theirs can be
   used with any S_N quadrature set, however, like the standard diamond scheme,
   *"our scheme is only first-order accurate."* This is a direct, primary-source
   statement on the BMC-43-vs-R&L-16 axis — see the one-line answer.
9. ⚠ **Notation:** τ is printed **without an m subscript** in both (16b) and
   (A2), though the right-hand side is manifestly m-dependent. Read `τ_m`
   everywhere; the paper never uses a level-uniform τ.
10. ⚠ **M&M's cylinder cosine convention differs from BMC's.** PDF p. 17 /
    pr. 630: μ = sin θ cos φ is the **radial** cosine and φ the azimuth; the
    ordinates lie on ξ-levels of common polar angle θ. This maps to ORPHEUS as
    M&M μ → `mu_x`, M&M φ → the azimuthal march variable ω. **Do not** carry
    M&M's letters into a BMC-notation context: BMC's radial cosine is η and its
    azimuth ω.

---

## Errata found (L-010 discipline: sidecar-vs-scan, then scan-vs-algebra)

| # | Where | Issue | Class |
|---|---|---|---|
| E1 | Eq. (1), PDF p. 4 / pr. 617 | The two ψ half-indices are **crossed**: as printed, `α_{m+1/2}` multiplies `ψ_{m−1/2}` and `α_{m−1/2}` multiplies `ψ_{m+1/2}`. Their own Eq. (4) — the same term with the ansatz substituted — is **uncrossed** (`α_{m+1/2}` with `μ_{m+1/2}`), as is Eq. (6a). The uncrossed form is the correct one: only it telescopes under the zeroth moment, which is what makes Eq. (5) the exact continuity equation. | **PRINT typo** (sidecar faithful; scan verified) |
| E2 | Eqs. (12) & (13), PDF p. 11 / pr. 624 | Sidecar OCRs the derivative as `dφ/dr`; the rendered page reads as `dψ/dr`, i.e. `dψ_s/dr`. `[M]` **the algebra is decisive**: with `dψ_s/dr`, Eq. (11)+(8) satisfies (12) exactly and (11)+(9) reproduces the `Q(1+2β)` of (13) to O(β); with `dφ/dr`, (11)+(8) yields 0 ≠ Q. Read `−dψ_s/dr`. | **OCR slip** (scan correct) |
| E3 | PDF p. 17 / pr. 630 | "μ = sin θ cos φ and **ξ = sin θ**". ξ is conventionally the *axial* cosine (cos θ); as printed it is the level's polar **sine**. Inert for (A2)–(A4) (which only use sin θ_ℓ and φ) and self-consistent as a level *label*, but do not import "ξ" from this paper into a Bell–Glasstone/BMC context. | **PRINT anomaly** (scan verified) |
| E4 | Summary, PDF p. 17 / pr. 630 | "…as long as the starting flux is not seriously underestimated relative to the **starting** fluxes." Should read "…relative to the **weighted** fluxes" (cf. the identical construction on pr. 629). Harmless. | **PRINT typo** |
| E5 | In-text ref. 2 vs the reference list | The text repeatedly credits "Miller and Alcouffe" for the origin-flux procedure, but reference 2 is **D. J. Dudziak, R. D. O'Dell, and R. E. Alcouffe**, "Transport and Reactor Theory," LANL report **LA-7911-PR** (July 1979) — no Miller in the author list. | **Attribution mismatch** |
| E6 | Cover page, PDF p. 1 | Sidecar OCRs the first author as "J. E. **Horel**" on the T&F metadata page. The article title page (PDF p. 2) prints **Morel**. | **OCR slip** (cosmetic) |

⭐ **Lead worth chasing (L-012 follow-up):** E5's **LA-7911-PR** (Dudziak,
O'Dell & Alcouffe, July 1979) is a strong candidate for the unresolved
"Alcouffe & O'Dell" primary behind Hébert's curvilinear cell-edge construction
(Hébert ref. 36), which four databases could not resolve and which remains
unread. It is a LANL **quarterly progress report** of Group T-1 — search OSTI,
not CrossRef. Status this session: the OSTI Python client hit
`ConnectionResetError` twice (network, not a miss); a web search confirms the
report exists and that the `LA-####-PR` T-1 series is hosted on OSTI in
full text (e.g. `osti.gov/servlets/purl/6513706`, `.../5521649`).
⛔ **UNVERIFIED and UNREAD** — I did not resolve the LA-7911-PR record itself,
and the search engine visibly conflated several quarters of the series (one hit
described an FY-81 quarter under a July-1979 label). Treat this as a lead only;
resolve the exact record before citing it. Note also that the *procedure* M&M
call "Miller–Alcouffe" is attributed in their reference list to a report with no
Miller author — so the true origin of the name is itself open.

---

## The one-line answer

**YES — the primary source resolves it, and it resolves it in favour of
BMC Eq. 43 (barycentric-in-μ), explicitly and by design.**

M&M pose τ as the barycentric coordinate of the node between its own cell's edge
cosines in *both* geometries (Eq. 16b sphere, Eq. A2 cylinder), never bound or
clamp it, and state in their own introduction — PDF p. 4 / pr. 617, verified on
the scan — that unlike the Reed–Lathrop scheme theirs works with **any** S_N
quadrature set but is, like plain diamond, **only first-order accurate**. The
BMC-43-vs-R&L-16 tension is therefore not an unresolved conflict to be
adjudicated: it is the **deliberate trade the primary source makes and
announces**. R&L buy second-order angular accuracy by *moving the ordinates*
(their Eq. 13b makes the quadrature an output, which M&M reject at pr. 616 on
the grounds that the resulting sets cannot integrate polynomials of degree > 3);
M&M buy diffusion-limit consistency (β = 0) on an arbitrary quadrature and pay
for it in order. `[M]` the price is visible in M&M's own default configuration —
on plain Gauss spheres τ₁ saturates at ≈ 0.3897 and |τ−½|/w diverges as O(N) at
the end ordinates, exactly R&L's degradation mechanism — and their cylinder
formula (A2)–(A4) on an equal-weight level yields τ₁ → ¼ / τ_M → ¾ with
(τ₁−½)/w = −M/4, which is precisely ORPHEUS's measured behaviour. **Nothing in
the primary licenses clamping τ toward ½; doing so would trade away the β = 0
property that is the paper's entire result, in exchange for an accuracy order
the paper never claims.**
