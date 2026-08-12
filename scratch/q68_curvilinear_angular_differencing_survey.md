# Q68 — Curvilinear SN angular differencing beyond the weighted diamond

**Status:** IN PROGRESS (written incrementally — a partial file is intentional).
**Author:** literature-researcher, 2026-08-12.

**Brief:** Q1 taxonomy of curvilinear angular-differencing schemes; Q2 what
production codes ship; Q3 is the (eta,phi) accuracy floor a named phenomenon;
Q4 does anyone measure diffusion-limit-consistency vs raw accuracy; Q5
post-2010 work.

**Verification legend**

- `[SCAN p.N]` — verified against the rendered PDF page N (1-based PDF page).
- `[OCR p.N]` — taken from the Mistral sidecar only, NOT yet page-verified.
- `[DB]` — bibliographic fact resolved against CrossRef/OpenAlex/OSTI.
- LOCAL / NOT LOCAL — presence in `scratch/literature/`.

---

## HEADLINE — the single most important source, and it was already local

**K. D. Lathrop, "A Comparison of Angular Difference Schemes for One-Dimensional
Spherical Geometry S_N Equations," Nucl. Sci. Eng. 134(3):239–264 (2000),
DOI 10.13182/NSE00-A2114. LOCAL** —
`scratch/literature/Lathropp(2000)A Comparison of Angular Difference Schemes for One-Dimensional Spherical Geometry SN Equations.pdf`
(sidecar mapping: printed p. = PDF p. + 237).

This paper *is* the Q1 taxonomy, and it answers Q4 with numbers. It compares
**five** angular-differencing schemes (plus a sixth, step, discussed and
dismissed) for the 1-D sphere, derived with **no spatial differencing at all**
— the angular discretization produces coupled ODEs in `r`, integrated with a
Mathematica ODE solver to 6 significant figures, and compared against a
**closed-form two-region purely-absorbing analytic solution** (his Eqs. 3, 5a–5d).
That construction isolates angular error exactly the way ORPHEUS's #235 floor
question needs.

⭐ **Lathrop's Eq. (30) + the sentence following it is the published statement of
the trade ORPHEUS measured.** `[SCAN p.8]` For a weighted diamond with node
offset `δ` (his Eq. 23; `δ = 2τ − 1`), the third Taylor coefficient of the
angular-derivative approximation is

```
 −(1/4) Δμ [ 2 δ (1 − μ_m²) + (1 − δ²) Δμ μ_m ]                (Eq. 30)
```

"…which is still of `O(δΔμ + Δμ²)`. The next term is `O(δΔμ² + Δμ³)` and so on,
**so that only with `μ_m = μ̄` (`δ = 0`) is the truncation order `O(Δμ²)`**."

⚠ The sidecar renders the prefactor as `−1/3`; the **printed page says `−1/4`**
`[SCAN p.8]`. (Same slip already recorded in lessons L-010.)

⟹ Any `τ ≠ ½` — including the Morel–Montry `τ` and the Reed–Lathrop `τ` —
**demotes the angular truncation order from second to first**, by the primary's
own analysis. Lathrop restates it twice more in Sec. IV `[SCAN p.11]`:
*"even if the second-order weighted diamond approximation corresponding to this
δ is used in the solution, the equation truncation error is not `O(Δμ²)`"*, and
*"Even if Reed's 'optimum' weighted diamond difference relations are used, the
truncation errors are `O(δΔμ + Δμ²)`, as shown earlier."* He also notes the
partial escape: *"δΔμ may be as small as `Δμ²`"* — true for Gauss–Legendre away
from the endpoints (his §III.B), which is exactly where it is NOT true for
ORPHEUS's grazing directions.

This is the mechanism behind ORPHEUS's `[M]` 2.5–3.0× — see §Q4.

---

## Q1 — The taxonomy of curvilinear angular-differencing schemes

### Q1.0 Lathrop 2000's five schemes (the spine of the taxonomy)

All notation is Lathrop's. Map to ORPHEUS: his `μ` = the marching cosine
(sphere: the streaming cosine; ORPHEUS cylinder: the radial cosine `η`), his
`Δμ_m` = the angular cell width, his `μ̄_m = (μ_{m+1/2}+μ_{m−1/2})/2` = the cell
midpoint, his `α_{m±1/2}` = ORPHEUS's α (his recursion Eq. 19 is
`α_{m+1/2} − α_{m−1/2} = −2 μ_m Δμ`, i.e. Lathrop–Carlson with `w_m = 2Δμ_m`
absorbed — check the factor-of-2 against ORPHEUS's `w`), and

```
    τ_ORPHEUS = (1 + δ_Lathrop)/2 ,   δ = 2τ − 1 ,   |δ| ≤ 1  ⟺  τ ∈ [0,1]
```

`[SCAN p.8]` — his Eq. (24) is literally ORPHEUS's closure with that
substitution. He writes the admissibility exactly as ORPHEUS's clamp question:
*"where `|δ| ≤ 1` as long as `μ_m` is in `Δμ`."* — i.e. the bound is
**node-in-cell**, nothing more; NO source in this survey clamps `τ` tighter
(the standing denominator is now **6 of 6 primaries**, see §Q4.3).

| # | Scheme | Defining equation(s) | Unknowns per cell | Angular order (Lathrop's own) | Cost |
|---|---|---|---|---|---|
| 0 | **Step / constant discontinuous** | Eq. (7) with `ψ_{m+1/2} = ψ_m`, expanded Eq. (8) | 1 | `O(Δμ)` | 1 ODE/cell |
| 1 | **Diamond** (plain, `δ=0`) | Eq. (18) + `ψ_{m+1/2} = 2ψ_m − ψ_{m−1/2}`; γ Eq. (21)/(22) | 1 | `O(Δμ²)` **iff `μ_m = μ̄`** | 1 ODE/cell |
| 2 | **Weighted diamond** (Reed–Lathrop "optimum") | Eq. (24) closure + Eq. (19) α-recursion + Eq. (28) node quadratic | 1 (+ the ordinates are OUTPUTS) | `O(δΔμ + Δμ²)` — 1st order unless `δ→0` | 1 ODE/cell |
| 3 | **Linear continuous** (the ORIGINAL S_N of Carlson) | Eq. (9) Taylor segment + Eq. (11)/(12) | 1 | `O(Δμ²)` | 1 ODE/cell |
| 4 | **Linear discontinuous in angle** (Walters & Morel 1991) | Eqs. (12) + (33), slope closed by Eq. (32) | **2** (`ψ_m`, `ψ_{m+1/2}`) | `O(Δμ²)` pointwise; **4th-order** on absorption/leakage `[SCAN p.19–20]` | ≥2× arithmetic per Δμ |
| 5 | **Quadratic continuous** (NEW in this paper) | Eqs. (38) + (41) [or (43)]; quadratic Eq. (35), slopes Eqs. (36)/(37) | **2** | all terms perfect-derivative or `O(Δμ⁴)`; **4th-order** measured | ≥2× arithmetic per Δμ |

Cost, in the author's own words `[OCR p.27]`: *"For a given N, the higher-order
schemes require at least twice as much arithmetic for each Δμ."* He explicitly
declines to cost them once space is also discretized.

**Derivation machinery** (this is the reusable part). Every scheme except
diamond is derived by *assuming an angular shape inside the cell and taking
normalized moments* `M_j = ∫ μ^j (·) dμ / Δμ` (his Eq. 6) of the **conservation
form** of the transport equation. He is emphatic about why: multiplying the
zeroth-moment equation by `Δμ` and summing telescopes the angular-derivative
term away, giving exact particle balance — *"a property of the analytic
transport equation, which is very useful to preserve"* `[OCR p.5]`. Diamond is
the odd one out: *"The diamond difference equations … cannot be derived by the
moment process. Instead, they are a hybrid scheme."* `[OCR p.7]`

**Scheme 5's structural novelty, and it is the most transferable idea in the
paper.** The second equation is obtained NOT by taking the first moment of the
transport equation but by taking the **zeroth moment of the transport
equation's first angular derivative** (his Eq. 34). Because every term of Eq.
(34) is already an angular derivative, its zeroth moment is *just the transport
equation evaluated at `μ_{m+1/2}` minus at `μ_{m−1/2}`, divided by Δμ* — so an
anisotropic Legendre source needs **no moment integrals at all**, only
`S_{m+1/2} − S_{m−1/2}`. `[OCR p.9]` He compares the two routes (his Eq. 39 vs
Eq. 43) and finds they differ in exactly two terms, with numerics giving "an
edge" to the derivative route.

**The `μ = −1` seed is a first-class concern in all five.** Abstract `[OCR p.2]`:
*"the desirability is shown of using an initializing computation of the `μ = −1`
angular flux to correctly compute the central flux and of having a difference
approximation that ensures this central flux is the same for all directions."*
Step and linear-discontinuous **decouple** from the `μ=−1` equation (their
`ψ_{m−1/2}` coefficient is `1 − μ²_{m−1/2} = 0` at `μ_{1/2} = −1`), which
propagates one wrong centre value into every subsequent direction — at `N=2` a
**36 % centre-flux error** `[OCR p.5]`. Walters & Morel's fix, adopted by
Lathrop, is a **hybrid**: run a *continuous* scheme in the first angular
interval only, then switch. (Directly relevant to the ORPHEUS "march SEED"
resume point.)

⛔ **Geometry caveat, load-bearing:** Lathrop 2000 is **SPHERE ONLY**. `[M]`
`grep -ci cylind` on the sidecar returns **0**. Every order claim above is for
the 1-D sphere with a *single* angular variable. Lifting schemes 4/5 to the
cylinder is not done in this paper and, per §Q1.4, is not done anywhere I found.

### Q1.1 The field's own verdict on the taxonomy — and it is a NEGATIVE

**E. W. Larsen and J. E. Morel, "Advances in Discrete-Ordinates Methodology,"
Ch. 1 of *Nuclear Computational Science: A Century in Review* (Y. Azmy and
E. Sartori, eds.), Springer (2010). LOCAL** —
`scratch/literature/Nuclear Computational Science - A Century in Review.pdf`
(printed p. = PDF p. − 12). **§1.5.1 "Angular Derivatives"** is a two-page
review of *exactly* this question.

⚠ **Self-correction.** I first wrote here that the existing memory entry for this
book (`larsen_morel_2010_review_extraction`) "does not record it". **That is
false** — the entry does log `pp.36-37 flux dip 1.92-1.93, ~5 pp` and `2-D cyl =
ONE sentence`. What the entry does NOT carry is the section's **verdict prose**,
which is the load-bearing part: the three quotes below. So the correct
characterisation is *"indexed but not mined"*, not *"missed"*.

`[SCAN p.49]` (printed p. 37), quoting tightly:

- *"Very few practical improvements beyond the elimination of the flux dip have
  been made in S_N angular derivative treatments."*
- *"Discontinuous finite-element discretizations might have been expected to
  have had an impact, but this has not happened, partly because it is difficult
  to develop a discontinuous angular finite-element method that is compatible
  with the standard S_N method in multidimensional geometries."*
- *"one of the few linear-discontinuous S_N angular derivative treatments ever
  developed for the 1-D spherical geometry equation was found to be **less
  accurate than the weighted-diamond scheme** … for a series of test problems
  [60]."* Their diagnosis: *"the starting-direction flux is computed … with
  greater accuracy than the other directions; hence, significant accuracy is
  actually lost if the starting-direction flux plays no role in the angular
  derivative treatment."*
- *"However, superior accuracy relative to the weighted-diamond scheme was
  obtained by using a **quadratic-continuous approximation in the first angular
  cell and a linear-discontinuous approximation in the remaining angular
  cells** [60]."*
- Closing: *"All of these factors make it challenging to develop advanced S_N
  angular derivative treatments [100]."* ([100] = Lathrop 2000.)

Their Eq. **(1.93)** `[SCAN p.49]` is the Morel–Montry weighted diamond in
ORPHEUS's own variable:

```
 ψ_n = ((μ_n − μ_{n−1/2})/w_n) ψ_{n+1/2} + ((μ_{n+1/2} − μ_n)/w_n) ψ_{n−1/2}
   ⟹  τ_n = (μ_n − μ_{n−1/2}) / w_n      (barycentric coordinate of the node)
```

and they state `[SCAN p.49]` *"This scheme has been generalized to 2-D
cylindrical geometry [38]"* — ref [38] is Morel & Montry 1984, i.e. the cylinder
appendix ORPHEUS already reproduces. **No other cylinder angular-differencing
generalization is named anywhere in the chapter.**

They also give the three-part diagnosis of the flux dip `[SCAN p.49]`, of which
the third is ORPHEUS's exact issue: *"The diamond-in-angle equation is
inconsistent with the location of the quadrature cosines within each angular
bin. As a result, the diamond-in-angle scheme does not preserve solutions that
are linear in μ."* And the seed fix is credited to **[26] = Dudziak, O'Dell &
Alcouffe, LA-7911-PR (1979)** — *"discretize the slab-geometry form of the
starting-direction flux equation rather than the curvilinear-like form, and set
all the angular fluxes at the centre of the sphere equal to the
starting-direction flux value."* (That closes the standing memory lead: LA-7911-PR
is the starting-direction primary. NOT LOCAL.)

⟹ **Q1 answer in one line: the published taxonomy is Lathrop 2000's five
schemes; everything else in 40 years is either (a) a variant of the weighted
diamond, or (b) an angular DFEM attempt that the field's own review says did not
pay off.**

### Q1.2 The linear-discontinuous-in-angle primary

**W. F. Walters and J. E. Morel, "Investigation of Linear-Discontinuous Angular
Differencing for the 1-D Spherical-Geometry S_N Equations," Proc. ANS Topl.
Mtg. *Advances in Mathematics, Computations, and Reactor Physics*, Pittsburgh,
28 Apr – 2 May 1991, Vol. III, p. 13.2 3-1 (ANS 1991). NOT LOCAL.**

⚠ Citation conflict between two local sources, both authored by people who were
there: Lathrop 2000 ref. 3 gives *"Vol. III, p. 13.2 3-1"*; Larsen–Morel ref.
[60] gives *"Sec. 11.1"*. Same meeting, same title. Treat the volume/section
pointer as unresolved until the proceedings are in hand.

This is the ONLY dedicated linear-discontinuous-in-angle paper for curvilinear
S_N that either review names. Its content is recoverable from Lathrop 2000
§III.D (Eqs. 12 + 32 + 33) since Lathrop re-derived and re-tested it, so
acquiring the proceedings is **optional, not blocking**. Its two documented
properties: (i) 4th-order convergence on integral quantities (absorption,
leakage) `[SCAN p.19–20]`; (ii) it **decouples from the μ = −1 seed**, which is
why it loses to weighted diamond pointwise unless hybridised.

### Q1.3 Angular schemes NOT found for curvilinear geometry — with search scope

Each of these was searched for explicitly. A negative with its scope:

| Sought | Searched | Result |
|---|---|---|
| **Step / constant discontinuous in ANGLE** | Lathrop 2000 §III (it is his Eqs. 7–8) | EXISTS but is `O(Δμ)` and decouples from the seed; Lathrop tested a "hybrid step + diamond first direction" and reports it *"displayed first-order error behavior in calculating absorption"* `[OCR p.20]`. Dead end, and it is a *measured* dead end. |
| **Characteristic / exponential differencing in ANGLE** | corpus grep (`characteristic`, `exponential`) + OpenAlex + Scopus + OSTI | **Found nothing for the angular variable.** Every hit (Walters–Wareing 1996 nonlinear characteristic, Askew 1972, Halsall 1980 CACTUS) is *spatial*. See §Q1.5 for the one adjacent case (Morel's Fokker–Planck angular scheme). |
| **DG / higher-order FEM in ANGLE for curvilinear S_N** | corpus + OpenAlex + Scopus | Only Walters–Morel 1991 (LD) and Lathrop 2000 (quadratic continuous). Larsen–Morel's explicit statement that angular DFEM "has not happened" is the authoritative denominator. |
| **Angular "nodal" schemes** | OpenAlex/Scopus/OSTI | The one hit — Wu & Roberts (1999) NSE 133 — is a **spatial** nodal method (see §Q5.1). |
| **2-D (η,φ) treatment of the cylindrical redistribution** | see §Q1.4 | **Nothing found.** |

### Q1.4 ⭐ The (η,φ) two-dimensional treatment

⛔ **REFUTED 2026-08-12, by me, within the same session.** I first wrote this
section as "searched hard, found NOTHING" on the strength of the corpus grep +
OpenAlex. That was wrong, and the reason it was wrong is instructive: **OpenAlex
full-text search is poisoned by `SN` = supernova** (three query formulations
returned WMAP, Type Ia spectropolarimetry and MHD papers). Switching to
CrossRef `journal_search` scoped to the NSE ISSN found the counterexample on the
first query. The original negative and its scope are kept below because the
*scope* is still the useful part; the *conclusion* is replaced by §Q1.4a.

**Original (refuted) scope of the negative:**

- Local corpus: greps for `two.dimensional angular`, `azimuthal derivative`,
  `phi derivative`, `both angular`, `omega derivative` across all 67 sidecars.
- OpenAlex: 3 query formulations (above) + `cylindrical geometry angular
  derivative two angular variables`.
- Scopus/CrossRef/OSTI/INIS: see §Q5.
- The two authoritative reviews (Larsen–Morel 2010 §1.5.1; Lathrop 2000) name
  no such scheme.

**What the primaries actually do instead**, and this is the important structural
point: Morel & Montry 1984's cylinder appendix (A1)–(A4) — already extracted by
ORPHEUS — handles the two angular variables by a **tensor-like split**, not a
2-D closure: partition in the AZIMUTH by cumulative weight, take the barycentre
in the RADIAL cosine. Hébert §3.9.3 does the same. So the "1-D η thread" is not
an ORPHEUS simplification of a richer published scheme — **it is what the
literature does.** The (η,φ) coupling is carried by the *α/W* azimuthal
streaming coefficients, not by a 2-D angular difference stencil.

⟹ *(original conclusion, now superseded)* "If ORPHEUS wants a genuinely 2-D
angular closure on the cylinder, there is no published scheme to port."

### Q1.4a ⭐⭐ There IS a published 2-D (μ,φ) angular closure — and it is NEW to this corpus

**F. Chaland and G. Samba (CEA/DAM/DIF), "Discrete Ordinates Method for the
Transport Equation Preserving One-Dimensional Spherical Symmetry in
Two-Dimensional Cylindrical Geometry," Nucl. Sci. Eng. 182(4):417–434 (2016),
DOI 10.13182/NSE15-38.** `[DB]` CrossRef-verified (authors, vol, pages, date).
**NOW LOCAL** — I extracted it from the user's own NSE archive
(`/Users/rodrigo/Downloads/NSE/Vol_182(4)_Nuclear Science and Engineering.zip`)
into `scratch/literature/Chaland-Samba(2016)…pdf` and OCR'd it (19 pp, $0.02).
Sidecar at `scratch/literature_ocr/Chaland-Samba(2016)….md`; printed p. = PDF p. + 415.

This is the only paper I found that discretizes **both** angular derivatives of a
cylindrical transport equation with an explicit difference stencil in each.

**The equation** `[SCAN p.6]`, their Eq. (6) — conservative form, with `μ`
measured against the vector from a *fixed point O on the z-axis* to the current
position (their "mixed coordinates", Appendix A), NOT against the z-axis:

```
 ∂_r[(μ r/R + √(1−μ²) cos φ · z/R) r u]  +  ∂_z[(μ z/R − √(1−μ²) cos φ · r/R) r u]
   + (1/R) ∂_μ[(1−μ²) r u]                         ← FIRST angular derivative
   − √(1−μ²) (z/R)(1/r) ∂_φ[sin φ · r u]           ← SECOND angular derivative
   + Σ r u  =  f r
```

**The angular closure** `[SCAN p.6]` — this is the part ORPHEUS wants:

| element | their choice |
|---|---|
| angular mesh | a **Cartesian product `(μ,φ)` mesh**, explicitly *"instead of using the standard nonconformal equi-solid-angle mesh"*. `μ_l = −1 + 2(l−1)/n`; `φ_m = π − 2π(m−1)/(n+2)`. Both equispaced. |
| closure, Eq. (8) | **plain diamond in BOTH variables** — `ū_{m+½,l+½} = (ū_{m+1,l+½} + ū_{m,l+½})/2 = (ū_{m+½,l+1} + ū_{m+½,l})/2`, i.e. `τ ≡ ½` twice over |
| μ-coefficients, Eq. (9) | `β_{l+1} − β_l = −2 μ̄_{l+½} Δμ`, `β_1 = 0`, **"note that `β_l = 1 − μ_l²`"** |
| φ-coefficients, Eq. (10) | `α_{m+1} − α_m = −cos φ̄_{m+½} Δφ`, `α_1 = 0`, `α_{m+1} − α_m ≈ sin φ_{m+1} − sin φ_m` |
| sweep order, Eq. (7) | from the angular characteristics: `dμ/ds = (1−μ²)/R > 0`, `dφ/ds = −√(1−μ²) sin φ (z/(Rr)) < 0` (for `z ≥ 0`) ⟹ **increasing `l` (μ: −1→+1), then increasing `m` (φ: π→0)** — a genuine 2-D angular sweep |
| spatial | finite volume, area-weighted (Wilkins-style); GMRES, **no sweep** *"since the characteristics are no longer straight lines"* |

⭐ **The `β_l = 1 − μ_l²` parenthetical is the cleanest published confirmation of
a fact ORPHEUS has been circling**: the conservation recursion lands *exactly* on
the geometric value `1 − μ²` **precisely when the node is the midpoint of an
equispaced-μ cell** (telescoping `Σ 2μ̄_j Δμ_j = Σ(μ²_{j+1} − μ²_j)`). With any
other node placement it does not, and the α/β sequence carries a defect — which
is Lathrop's `β_{m+1/2}` of his Eq. (26).

**Claimed order and cost.** Abstract: *"it is consistent of order 1"* — first
order, and the authors are candid that the limitation is the *spatial* scheme,
not the angular one (§V.C uses MMS with the angular mesh refined until *"the
error due to the angular discretization is negligible"*). Cost: they abandon
sweeping entirely for GMRES.

⛔ **Two disqualifiers before anyone tries to port this to ORPHEUS's 1-D cylinder:**
1. **`[SCAN p.16]` their own words: *"our scheme does not preserve the diffusion
   limit."*** They fix it in Appendix B by swapping in Stone finite elements
   (refs 14/15) — *"strict 1-D symmetry and order 2 is obtained"* — and defer the
   diffusion-limit proof to a follow-up paper.
2. The whole construction is **motivated by ray effects in 2-D ICF hydro**, not
   by 1-D accuracy, and the mixed-coordinate frame is *"no longer invariant under
   translation along the z-axis"*, so *"ray effects may be present for
   calculations far from the 1-D spherical symmetry, in which case the standard
   form of the equation should be retained."*

⟹ **Corrected Q1.4 answer.** A 2-D `(μ,φ)` product-mesh angular closure with
plain diamond in both variables and two independent conservation recursions
**does exist in print, is recent, and is in NSE.** What does NOT exist is such a
closure for the *standard* 1-D/2-D cylindrical `(η,ξ)` frame with a stated
angular convergence order and a diffusion limit. Chaland–Samba is the structural
template; it is not a drop-in.

**Follow-up worth pulling** (NOT LOCAL, not yet resolved): the "following paper"
they promise on the diffusion limit of the Stone-FE variant. Also their ref. 13
(the "discrete element method" they say their frame-adaptation resembles) and
refs 14/15 (Stone finite elements) — I did not resolve these; the reference list
is on PDF p. 17.

**Original scope note retained:** confidence high for "not in the two reviews"
(Larsen–Morel 2010 and Lathrop 2000 genuinely do not name any 2-D angular
closure); the M&C/PHYSOR proceedings 1991–2010 remain unexhausted.

---

## Q4 — The diffusion-limit vs raw-accuracy trade: YES, it is measured, and by Lathrop

**Answer: the trade is real, published, quantified, and the sign and magnitude
of ORPHEUS's `[M]` 2.5–3.0× match the published data. This is not a signal that
the ORPHEUS implementation is wrong.**

Three independent legs.

### Q4.1 The MECHANISM is stated by the primary

Lathrop 2000 Eq. (30) `[SCAN p.8]`, already quoted in the headline: with node
offset `δ = 2τ − 1`, the third Taylor coefficient is `O(δΔμ + Δμ²)`, and
**"only with `μ_m = μ̄` (`δ = 0`) is the truncation order `O(Δμ²)`."** Sec. IV
repeats it for the specific case ORPHEUS cares about `[SCAN p.11]`: *"Even if
Reed's 'optimum' weighted diamond difference relations are used, the truncation
errors are `O(δΔμ + Δμ²)`, as shown earlier."*

`[SCAN p.11]` also pins the two-unknown schemes' orders, which §Q1.0's table
depends on: for linear discontinuous the three terms carry `O(Δμ⁴)`, `O(Δμ²)`,
`O(Δμ⁴)` errors (spatial derivative / angular derivative / cross-section), so
*"the overall accuracy of the linear discontinuous approximation is `O(Δμ²)`"*;
for quadratic continuous *"all the terms are either perfect derivatives or give
errors of `O(Δμ⁴)`."* And his honest caveat: *"While such error analysis is
useful, it does not give the whole picture. The linear discontinuous
approximation gives better results than predicted for integral quantities such
as the absorption rate."*

⟹ **Plain diamond with a midpoint node is second-order in angle. Every weighted
diamond — Reed–Lathrop's and Morel–Montry's alike — is first-order in angle,
with the first-order coefficient proportional to `δ`.** The diffusion-limit
property is bought with exactly that order.

Corroborated from the other side by Lathrop's own §V `[OCR p.12–13]`: plain
diamond with `μ_m = μ̄` *"would not have `c = 2/3` and hence would not give the
correct diffusion limit"* (his `c ≡ Σ Δμ_m μ_m²`, Eq. 54; on an equal-interval
mesh `c = 2/3 − Δμ²/6`, his Eq. 63). **So the two properties are in direct
opposition on a fixed quadrature: `δ = 0` buys `O(Δμ²)` and loses `c = 2/3`;
`δ ≠ 0` buys `c = 2/3` and loses the order.** That is the trade, stated by name.

### Q4.2 The trade is MEASURED, direction by direction

Lathrop's Table II (maximum relative angular flux error, %, test problem 1;
`[SCAN p.14]`, values reproduced exactly from the rendered page). Ratio
`WD / diamond` computed by me `[M]`:

| N | μ | diamond | weighted diamond | WD/DD |
|---|---|---|---|---|
| 8 | +0.125 | 5.58 | 7.09 | **1.27** |
| 8 | +0.375 | 0.56 | 1.39 | **2.48** |
| 16 | +0.0625 | 3.51 | 4.00 | **1.14** |
| 16 | +0.1875 | 0.49 | 1.05 | **2.14** |
| 16 | +0.3125 | 0.14 | 0.55 | **3.93** |
| 16 | +0.4375 | 0.14 | 0.26 | **1.86** |
| 16 | +0.5625 | 0.14 | 0.15 | **1.07** |

⭐ **Weighted diamond is WORSE than plain diamond by up to 3.93×, and the number
of directions where it loses GROWS with N: 2 of 8 at S8, 5 of 16 at S16.** The
losing band is the **small-positive (outgoing near-grazing) cosines** — exactly
where ORPHEUS's cylinder azimuthal thread is weakest. ORPHEUS's `[M]` 2.5–3.0×
sits inside this published band.

⚠ **Read the table, not the abstract.** Lathrop's own summary sentence is
*"the weighted diamond difference approximation has smaller maximum and average
relative flux errors than the diamond or the linear continuous"* `[OCR p.18]`,
and on the *average* error and on most directions that is true. The
direction-resolved maximum error inverts it in the outgoing grazing band, and the
inversion strengthens under angular refinement. Both statements come from the
same table.

### Q4.3 The clamp question, answered again — 6 of 6 primaries, none clamps τ

Lathrop 2000 states the admissibility bound exactly as node-in-cell and nothing
more: *"where `|δ| ≤ 1` as long as `μ_m` is in `Δμ`"* `[SCAN p.8]`, i.e.
`τ ∈ [0,1]`. Standing denominator for "does any primary clamp τ away from the
endpoints":

| # | primary | clamps τ? |
|---|---|---|
| 1 | Reed & Lathrop 1970 | no (`τ ∈ [0.44,0.56]` is an *observed* S8 range, not a rule) |
| 2 | Morel & Montry 1984 | no ("any value between zero and one") |
| 3 | Bailey, Morel & Chang 2010 | no |
| 4 | Hébert 2009 | no (τ ≡ ½ throughout) |
| 5 | Larsen & Morel 2010 §1.5.1 | no (Eq. 1.93 is the bare barycentric formula) |
| 6 | **Lathrop 2000** | **no** — `|δ| ≤ 1` ⟺ `τ ∈ [0,1]`, node-in-cell only |

### Q4.4 ⛔ What Q4 does NOT license

- Lathrop's WD is the **Reed–Lathrop** variety, where the ordinates are OUTPUTS
  of Eq. (28), on an **equal-interval** μ mesh. ORPHEUS ships the
  **Morel–Montry** variety (ordinates are INPUTS from the quadrature). They are
  different schemes (lessons L-013). What transfers is the *mechanism*
  (`δ ≠ 0` ⟹ first order) and the *qualitative* direction-resolved inversion —
  NOT the specific ratios.
- Lathrop's problems are **purely absorbing spheres with isotropic sources**,
  deliberately: *"purely absorbing systems with the most rapidly varying angular
  flux provide the greatest challenge. As scattering is introduced … the
  numerical difference between results for the various approximations can be
  expected to decrease."* `[OCR p.15]` ORPHEUS's fixture is an **anisotropic
  cylinder MMS**, which is a different corner.
- ⚠ **Nobody in this survey measures the trade in the regime that would settle
  ORPHEUS's choice** — i.e. a diffusive cylinder where the WD's `c = 2/3` should
  pay off — against a non-diffusive one, on the *same* code. Lathrop's diffusion
  analysis is purely analytic (§V); his numerics are purely absorbing. **That
  specific measurement is a gap.**

---

## Q2 — What production codes ship

⛔ **Constraint declared up front.** My operating rules forbid referencing
export-controlled codes — names, manuals, or equation numbers — in any output,
and several of the codes named in the brief are restricted-distribution
(RSICC-controlled or vendor-proprietary). **I will not characterize their
internals**, even where I could infer them. What follows is built only from the
**openly published record**. Where a restricted code is the only possible source
for an answer, I say the answer is unavailable to me rather than guessing. If
the user needs those specific codes characterized, that has to come from their
own licensed documentation, not from me.

### Q2.1 The openly-documented production answer: plain diamond, `τ ≡ ½`

**A. Hébert, *Applied Reactor Physics*, Presses Internationales Polytechnique
(2009), Ch. 3. LOCAL.** Hébert's Ch. 3 is the published method description of an
openly-distributed lattice-code lineage, and it is unambiguous. For the **1-D
cylinder** (§3.9.3), Eq. **(3.406)** `[OCR]`:

```
 φ_{p,q,i} = ½(φ_{p,q,i−1/2} + φ_{p,q,i+1/2}) = ½(φ_{p,q−1/2,i} + φ_{p,q+1/2,i})
```

— plain diamond in **space** and in **angle**, `τ ≡ ½`, on the equal-weight
equispaced-ω quadrature of his Eqs. (3.369)–(3.376), with the α recursion
(3.399) `α_{p,q+1/2} = α_{p,q−1/2} + W_{p,q} μ_{p,q}` and the `ω = π`
starting-direction seed (3.407)–(3.410).

⭐ Note the *identical structural form* to Chaland–Samba's Eq. (8) `[SCAN p.6]`:
a single cell-average tied to BOTH a spatial mean and an angular mean by
midpoint averaging. Two independent code lineages, 7 years apart, one for
lattice physics and one for ICF radiative transfer, both choose **`τ = ½`**.

⟹ **The one production curvilinear-SN angular closure I can document openly uses
the scheme ORPHEUS measured as the BETTER one.** That is a meaningful
corroboration of the `[M]` 2.5–3.0×: ORPHEUS's "unprincipled" fallback is what a
mature production code ships.

### Q2.2 What the field's review says about production practice generally

Larsen & Morel 2010 §1.5.1 `[SCAN p.49]`, on the traditional
β-coefficients + diamond-in-angle + starting-direction treatment: *"which is
representative of the traditional treatment used in essentially all curvilinear
geometries."* And on the weighted-diamond upgrade: *"This scheme has been
generalized to 2-D cylindrical geometry."*

Lathrop 2000 `[OCR p.25]` gives the reason the *Reed–Lathrop* weighted diamond
did not become production practice, and it is not an accuracy reason:

> *"This scheme has not been widely used because it specifies a quadrature that
> may not integrate a high-enough-order polynomial to ensure a correct
> calculation of the zeroth moment of an anisotropic source, so balance would not
> be preserved."*

He then notes the same difficulty afflicts **all five** schemes, because the four
non-WD ones are most accurate at `μ_m = μ̄`, *"and this quadrature correctly
integrates only the zeroth and first powers of μ"* — and offers the fix, his
Eq. (68): redefine `P_k(μ_m)` as the **cell-averaged** Legendre polynomial
`∫ P_k dμ / Δμ = [(μP_k − P_{k−1})/((k+1)Δμ)]` evaluated at the cell edges. This
*"separate[s] the choice of quadratures for the streaming term and the source
term, allowing them to be optimized separately"*, and he credits the
interpolating-polynomial variant to **Morel 1989 NSE 101:72 (LOCAL)**.

⭐ **This is a concrete, cheap, actionable item for ORPHEUS**, independent of the
τ question: if ORPHEUS moves to a midpoint-node equispaced rule to recover
`O(Δμ²)`, Lathrop's Eq. (68) is the published remedy for the anisotropic-source
balance loss that move causes. It is the missing third piece that makes
"τ = ½ + equispaced" a *complete* scheme rather than a fallback.

### Q2.3 What I could NOT determine

- **Any restricted-distribution production code's angular differencing** —
  declined by constraint (see above), not by failure to find.
- Whether any production code documents the **(η,φ) floor** or a fix — searched
  the open literature (see Q3); **found nothing**.

---

## Q3 — Is the floor a KNOWN, NAMED phenomenon?

**Answer: NO for the azimuthal/(η,φ) accuracy floor specifically. There is no
name and no canonical citation.** But there are **three named neighbours**, and
one of them is much closer than it looks.

### Q3.1 Named phenomena that are NOT it

| name | canonical citation | what it actually is | why it is not the floor |
|---|---|---|---|
| **discrete-ordinates flux dip** | Morel & Montry 1984 TTSP 13(5):615 (LOCAL) | erroneous suppression of the flux at the **centre** of a diffusive sphere | a centre-point artefact from seed + node-inconsistency, *fixed* by the WD; ORPHEUS's floor persists with the WD |
| **ray effects** | Lathrop, "Ray Effects in Discrete Ordinates Equations," NSE 32:357 (1968) `[DB]`; Briggs–Miller–Lewis NSE 57:205 (1975); Miller–Reed NSE 62:391 (1977) | absence of direction-to-direction coupling in multi-D ⟹ streaked solutions | an *angular-quadrature* resolution artefact, not a differencing-order defect |
| **thick-diffusion-limit failure** | Larsen–Morel–Miller 1987 JCP 69 (LOCAL); Larsen–Morel 1989 JCP 83 (LOCAL) | a *spatial* discretization losing the diffusion limit at large `Σh` | spatial, not angular |

### Q3.2 ⭐ The closest published concept, and it IS a named mechanism

**J. Hu and Y. Y. Azmy** — two papers, both `[DB]` CrossRef-verified, **NOT LOCAL**:

1. "Asymptotic convergence of the angular discretization error in the scalar flux
   computed from the particle transport equation with the method of discrete
   ordinates," **Ann. Nucl. Energy 138:107199 (2020)**,
   DOI `10.1016/j.anucene.2019.107199`. **OPEN ACCESS** —
   `https://www.osti.gov/servlets/purl/1801122` (OSTI full text).
2. "On the Regularity Order of the Pointwise Uncollided Angular Flux and
   Asymptotic Convergence of the Discrete Ordinates Approximation of the Scalar
   Flux," **Nucl. Sci. Eng. 195:598–613 (2021)**,
   DOI `10.1080/00295639.2020.1860634`. Not OA; NSE vol. 195 should be in the
   user's archive.

Their named mechanism: **the "regularity order with respect to the azimuthal
angle" of the pointwise angular flux caps the achievable angular convergence
rate.** From the 2021 abstract `[DB via OpenAlex]`: a theory of the regularity
order in the **azimuthal** angle is derived from the integral transport equation;
it *"can be estimated for a given problem configuration"*; and it motivates a
**Modified Simpson's quadrature** that *"converges the uncollided scalar flux
faster than any of the traditional quadratures **by avoiding integration across
points of irregularity in the azimuthal angle**."* Different quadrature types are
measured converging at *different orders* on the same problem.

⚠ **Scope caveat, and it is the reason this is a lead and not an answer:** their
derivation is **2-D Cartesian**, and the mechanism is on the **quadrature**
(integration of a non-smooth integrand), not on the **angular difference
stencil**. So it is not a citation for "the cylindrical weighted-diamond
azimuthal floor". It *is* the right vocabulary — "azimuthal regularity order
limits the observed convergence rate" — and their diagnostic (estimate the
regularity order, then check whether the quadrature straddles an irregularity)
is directly applicable to ORPHEUS's azimuthal thread.

### Q3.3 Search scope for the Q3 negative

So that the next session does not repeat it. `[M]` all run 2026-08-12:

- **CrossRef `journal_search`**, 4 query formulations × 3 ISSNs (NSE `0029-5639`,
  ANE `0306-4549`, JCP `0021-9991`) = 12 searches: `azimuthal angular
  discretization error cylindrical transport`, `convergence order loss
  curvilinear discrete ordinates angular`, `flux dip discrete ordinates`,
  `accuracy limitation angular redistribution cylindrical Sn`. Plus 3 further
  formulations × 5 ISSNs (adding TTSP `0041-1450`, JQSRT `0022-4073`).
- **OpenAlex**, 7 free-text formulations. ⚠ **Low yield — `SN` collides with
  "supernova" and swamps the results.** Do not rely on OpenAlex full-text for
  `S_N` topics; use ISSN-scoped CrossRef.
- **Scopus** — unavailable this session (`403 Forbidden` on every query; API key
  or network). **NOT searched.** This is the one real hole in the denominator.
- **Local corpus** — 67 sidecars grepped for `azimuthal derivative` (0 hits),
  `two-dimensional angular` (1 hit, LA-3251), `exponential angular` (0),
  `characteristic angular` (0), `angular nodal` (0).
- **Zotero** — not consulted (see §Tooling note at the end).

⟹ **Reporting recommendation: do not tell a reader the floor is a known named
phenomenon.** Describe it mechanically instead — *"the angular truncation error
of a weighted-diamond closure is `O(δΔμ + Δμ²)` (Lathrop 2000, Eq. 30), so it is
first-order in angle whenever the ordinate is not the cell midpoint"* — which is
citable, checkable, and true. That is the L-012 "restate as a MECHANISM, not a
survey" move.

---

## Q5 — Post-2010

**Answer: there is post-2010 work on curvilinear S_N, but almost none of it is on
the ANGULAR difference scheme. The one genuine exception is Chaland & Samba 2016
(§Q1.4a). BMC 2010 really is still the last word on angular-differencing
diffusion-limit analysis.**

Method: the reliable net was the **citation graph** (Semantic Scholar
`get_citations` on Morel–Montry 1984 and Lathrop 2000) plus ISSN-scoped CrossRef.
`[M]` Morel–Montry has ~30 indexed citers; Lathrop 2000 has **4**.

### Q5.1 Post-2010 hits, classified

| ref | what it actually is | angular scheme? |
|---|---|---|
| **Chaland & Samba (2016) NSE 182:417–434**, `10.13182/NSE15-38`. **LOCAL** (I added it) | 2-D `(μ,φ)` product-mesh angular closure, plain diamond in both, mixed-coordinate cylindrical frame | ⭐ **YES** — see §Q1.4a. The one real hit. |
| **Sun, Jiang & Xu (2019)**, "Multiscale Radiative Transfer in Cylindrical Coordinates," *Commun. Appl. Math. Comput.* 1:117–139, `10.1007/s42967-019-0007-x`. NOT LOCAL | unified gas-kinetic scheme in cylindrical coords; cites Morel–Montry | ⚠ **UNRESOLVED** — abstract not retrievable from OpenAlex/CrossRef. Cites M&M, so it must handle the angular derivative. **Worth one follow-up pull.** |
| **Liu, Tan & Wei (2026)**, "Taylor basis DG scheme for S_N transport equations on non-orthogonal grids in r-z coordinates," *Ann. Nucl. Energy* 233:112268, `10.1016/j.anucene.2026.112268`. NOT LOCAL. Preprint twin: `10.2139/ssrn.5447817` | DG on r-z grids — brand new (Aug 2026) | likely **SPATIAL** DG; r-z S_N carries the angular derivative so its treatment is stated somewhere. **UNRESOLVED.** |
| **Dean Wang (2019) NSE 193:1339–1354**, `10.1080/00295639.2019.1638660` | *"The Asymptotic Diffusion Limit of Numerical Schemes for the S_N Transport Equation"* — result: the mesh size needed to attain the diffusion limit is `ε^{1/k} h` with `k` the order of **spatial** accuracy `[DB abstract]` | ⛔ **NO — spatial.** Despite the title's resemblance to BMC 2010, this is the *spatial* diffusion limit. Do not confuse the two. |
| **Hu & Azmy (2020 ANE / 2021 NSE)** | azimuthal **regularity order** caps angular convergence; Modified Simpson's quadrature | ⛔ quadrature, 2-D Cartesian — see §Q3.2. Closest concept, wrong object. |
| **Machorro (2007) JCP 223:67–81**, `10.1016/j.jcp.2006.08.020`. NOT LOCAL | DG-FEM for the non-scattering 1-D **spherical** transport equation, various trial/test spaces, compared against diamond differencing and corner balancing | ⛔ **SPATIAL** DG. Included because the title reads angular and is not. |
| **Larsen (2003) TTSP 32:623–643**, `10.1081/TT-120025069` | *"Infinite-Medium Solutions of the Transport Equation, S_N Discretization Schemes, and the Diffusion Approximation"* | spatial/infinite-medium; not curvilinear angular |
| ANE 2023 wavelet angular discretization; NET 2023 hp-angular adaptivity; PNE 2024 goal-based angular adaptivity; NSE 2015 Ahrens *Lagrange Discrete Ordinates*; NSE 2017 DFEM-based quadratures; JCP 2023 FEM angular discretization on spherical geodesic grids | all **angular-QUADRATURE / angular-adaptivity** work in **Cartesian / multi-D** | ⛔ none treats a curvilinear angular *derivative* |

### Q5.2 The post-2010 negative, with scope

`[M]` 2026-08-12. **Searched and found nothing further:**
- Semantic Scholar citation graph of Morel–Montry 1984 (all ~30 citers listed and
  classified) and of Lathrop 2000 (all 4 citers: NSE 2004 singular solutions;
  NSE 2009 new boundary condition; *J. Fusion Energy* 2011 discontinuous 1-D
  spherical solvers w/ diffusion preconditioner; NSE 2012 multiproblem strategy —
  **none is a new angular scheme**).
- CrossRef ISSN-scoped searches over NSE / ANE / JCP / TTSP / JQSRT (see §Q3.3).
- **NOT searched:** Scopus (403 all session), INIS, OSTI, HAL, J-STAGE, Zenodo,
  and the M&C / PHYSOR / ANS-Transactions proceedings series. **The proceedings
  gap is the significant one** — Walters & Morel 1991, the single most important
  LD-in-angle primary, is itself a proceedings paper, which is direct evidence
  that this sub-field publishes there.

⟹ If the user wants the post-2010 denominator closed, the next move is
**proceedings**, not journals: M&C 2011/2013/2015/2017/2019/2021/2023,
PHYSOR 2012–2024, and ANS Transactions. I did not search these and they are
poorly indexed by the tools available here.

---

## Consolidated answer, and the three things I would act on

1. **⭐ Read Lathrop 2000 §III–V and Table II before anything else.** It is
   LOCAL, it is the taxonomy, and its Eq. (30) explains ORPHEUS's `[M]` 2.5–3.0×
   as a *published, expected* trade rather than a bug: any `τ ≠ ½` is
   **first-order in angle**; only the midpoint node is second-order. His own
   Table II shows weighted diamond losing to plain diamond by up to **3.93×** in
   the outgoing near-grazing directions, with the losing set growing under
   angular refinement. **ORPHEUS's measurement is consistent with the
   literature.**

2. **⭐ The `τ = ½` route has a known missing third piece, and Lathrop supplies
   it.** A midpoint-node equispaced rule integrates only `P_0` and `P_1`
   correctly, which breaks anisotropic-source balance. His Eq. (68) — replace
   `P_k(μ_m)` with the **cell-averaged** `∫P_k dμ/Δμ` — fixes it and decouples
   the streaming quadrature from the source quadrature. Credit for the
   interpolating-polynomial variant: **Morel 1989 NSE 101:72, LOCAL**. This makes
   "plain diamond + equispaced" a complete, defensible scheme rather than a
   fallback. It also happens to be **what the one openly-documented production
   code does** (Hébert Eq. 3.406).

3. **⭐ If ORPHEUS wants real order, the published upgrade path is
   two-unknowns-per-angular-cell**, and it is measured: linear discontinuous
   (Walters & Morel 1991) and quadratic continuous (Lathrop 2000 §III.E) both
   give **4th-order** absorption/leakage convergence at **≥2× arithmetic per
   Δμ**. Lathrop's quadratic continuous is the better of the two on pointwise
   flux error and does not decouple from the `μ=−1` seed. Its structural trick —
   take the zeroth moment of the transport equation's **first angular
   derivative** (his Eq. 34) instead of the first moment — makes the anisotropic
   source a bare difference `S_{m+1/2} − S_{m−1/2}` with no moment integrals.
   ⚠ But heed Larsen & Morel's warning `[SCAN p.49]`: the LD variant was
   *measured less accurate* than weighted diamond because it abandons the
   accurate starting-direction flux; the hybrid (continuous first cell, then LD)
   is what wins.

**And the one thing NOT to do:** do not treat the (η,φ) floor as a solved
problem with a scheme to port. Nothing in the literature gives a curvilinear
angular closure with both a stated angular convergence order and a diffusion
limit. Chaland–Samba 2016 has the 2-D angular structure but explicitly lacks the
diffusion limit.

---

## Bibliography with acquisition status

**LOCAL and read this session**
- Lathrop, NSE **134**(3):239–264 (2000), `10.13182/NSE00-A2114`. ⭐ the taxonomy.
- Larsen & Morel, "Advances in Discrete-Ordinates Methodology," Ch. 1 §1.5.1 of
  *Nuclear Computational Science: A Century in Review* (Azmy & Sartori, eds.),
  Springer. `[DB]` DOI **`10.1007/978-90-481-3411-3_1`**, CrossRef date
  **2009-12** (the book is usually cited as 2010). ⭐ the field's own verdict.
- Chaland & Samba, NSE **182**(4):417–434 (2016), `10.13182/NSE15-38`. ⭐ **newly
  added to `scratch/literature/` + OCR'd this session** from the user's NSE archive.
- Wu, Xie & Fischer, NSE **133**(3):350–357 (1999), `10.13182/NSE99-A2095`. **newly
  added + OCR'd.** ⛔ **SPATIAL** nodal Green's-function method; the angular
  treatment is the standard `β_{m+1} − β_m = −j W_m μ_m` recursion (their Eq. 15,
  `j = 1` cylinder / `j = 2` sphere) with the redistribution moved into the
  source. **Not an angular scheme** — it does bear on the *spatial* coarse-mesh
  question, which is a different ORPHEUS thread.
- Hébert, *Applied Reactor Physics* Ch. 3 (2009). Cylinder Eq. (3.406) = plain
  diamond in space AND angle.

**NOT LOCAL — ranked by value if the user wants them**
1. **Walters & Morel (1991)**, M&C Pittsburgh — the LD-in-angle primary.
   ⚠ two conflicting page pointers (Lathrop: Vol. III p. 13.2 3-1; Larsen–Morel:
   Sec. 11.1). **Low urgency** — Lathrop 2000 §III.D re-derives and re-tests it.
2. **Sun, Jiang & Xu (2019)** CAMC 1:117–139 — cylindrical UGKS, cites M&M,
   angular treatment unresolved.
3. **Liu, Tan & Wei (2026)** ANE 233:112268 — brand-new r-z DG; free SSRN
   preprint twin at `10.2139/ssrn.5447817`.
4. **Hu & Azmy (2021)** NSE **195**:598–613 — azimuthal regularity order.
   **Should be in the user's NSE archive** (`Vol_195(*)`). The 2020 ANE companion
   is **free** at `https://www.osti.gov/servlets/purl/1801122`.
5. **Dudziak, O'Dell & Alcouffe, LA-7911-PR (1979)** — the starting-direction-flux
   primary, now confirmed by Larsen–Morel ref. [26]. Closes a standing memory lead.
6. Machorro (2007) JCP 223:67–81 — spatial DG on the sphere; low priority.

**Bibliographic corrections `[DB]`**
- Bailey, Morel & Chang is **NSE 165:149–169 (2010-06), DOI `10.13182/NSE08-66`**
  — *not* `NSE08-64`. Worth checking `docs/refs.bib`.
- Lathrop 2000 is **NSE 134(3):239–264**, DOI `10.13182/NSE00-A2114`.

---

## Tooling notes for the next session

- ⛔ **Do NOT use OpenAlex free-text for `S_N` topics.** `SN` matches
  *supernova*; 7 formulations returned WMAP, Type Ia spectropolarimetry, MHD and
  sea-ice papers. **ISSN-scoped CrossRef `journal_search` is the right instrument**
  and found the key 2016 paper on the first try. This mistake cost me a false
  negative that I had to retract in §Q1.4.
- ⛔ **Scopus returned `403 Forbidden` on every query this session.** Not
  searched; this is the one acknowledged hole in every negative above.
- ⚠ Semantic Scholar rate-limits (`429`) after ~2 rapid calls; space them.
- ⭐ **The user's NSE archive is at `/Users/rodrigo/Downloads/NSE/`**, as
  `Vol_NNN(I)_Nuclear Science and Engineering.zip` bundles with per-article PDFs.
  `unzip -l` the bundles and grep the listing by title — this is how I found
  Chaland–Samba and Wu. **This is a much richer Tier-0 surface than
  `scratch/literature/` alone and it was not in my memory.**
- Zotero: **not consulted this session** — everything resolved from the local
  archive + CrossRef, so no annotations were checked. If τ-clamp provenance ever
  needs the user's own highlights, Zotero is still the place.
</content>
</invoke>
