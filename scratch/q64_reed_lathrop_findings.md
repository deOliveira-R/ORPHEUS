# Reed & Lathrop (1970) — the ANGULAR difference scheme

**Full citation.** Wm. H. Reed and K. D. Lathrop, "Truncation Error Analysis of
Finite Difference Approximations to the Transport Equation," *Nuclear Science
and Engineering* **41**(2), 237–248 (1970). DOI `10.13182/NSE70-A20710`.
Los Alamos Scientific Laboratory. Received 1970-02-06, revised 1970-03-30.

**Sources used.** PDF `scratch/literature/07-Truncation Error Analysis of Finite
Difference Approximations to the Transport Equation.pdf` (13 pages) + sidecar
`scratch/literature_ocr/<same stem>.md`. Page mapping **printed = PDF + 235**
(verified on the rendered headers: PDF p. 3 = printed 238 … PDF p. 13 = printed
248). Table I body recovered from the `.mocr.json` cache (sidecar carries only a
placeholder link — L-008).

**Scan verification.** Every equation quoted below was read on the rendered page
(PDF pp. 3–7 opened visually), not only in the sidecar. Two classes of defect
found and separated:

- ⚠ **OCR slips (sidecar wrong, print right).** The sidecar systematically
  misreads the half-integer subscript `1` as `3` on PDF p. 3: it renders Eq. (2)
  with `α'_{m+3/2}` and `N_{m−3/2}`, Eq. (4a) with `N_{i+3}`, and Eq. (4b) with
  `N_{m+3/2}`/`N_{m−3/2}`. **The print is correct** (`m±1/2`, `i+1`). Anyone
  grepping the sidecar for these equations must not transcribe the `3/2`s.
- ✅ **No published typo found.** Unlike BMC 2010 Eq. (52) (L-010), every
  load-bearing equation here is internally consistent and non-vacuous on the
  page. The sphere→cylinder analogue cross-read (Eqs. 13 vs 28) closes.

**Zotero.** Not consulted — the brief named a local file and the extraction is
complete from it; no annotations were checked.

---

## Q1 — What exactly is "Reed and Lathrop's weighted diamond scheme"?

### The defining equations

Printed p. 238 (PDF p. 3), sphere. The scheme has **two halves that share one
weight each**, and this pairing is the whole idea:

- **Where the node sits inside its cell** — Eqs. (3a)/(3b):

  `r_{i+1/2} = a_{i+1/2} r_{i+1} + (1 − a_{i+1/2}) r_i`
  `μ_m = τ_m μ_{m+1/2} + (1 − τ_m) μ_{m−1/2}`

  attributed in the text to **Grant, J. Comp. Phys. 2, 381 (1968)** (their ref. 6).

- **The closure that reconstructs the outgoing face value** — Eqs. (4a)/(4b),
  *with the same two weights*:

  `N = a_{i+1/2} N_{i+1} + (1 − a_{i+1/2}) N_i`
  `N = τ_m N_{m+1/2} + (1 − τ_m) N_{m−1/2}`

Eq. (4b) is **exactly ORPHEUS's angular closure**, and Eq. (3b) is exactly the
barycentric-coordinate definition of τ (BMC 2010 Eqs. 42/43, my predicate P2).
So the *skeleton* is identical to what ORPHEUS ships. **The letter τ is even the
same letter.** What differs is which equations are allowed to determine it — see
Q6.

**Range statement, verbatim (printed p. 238):** the parameters `a_{i+1/2}` and
`τ_m` "are restricted to the interval (0,1)" so that the point `(r_{i+1/2}, μ_m)`
"be contained within the mesh cell under consideration, which seems desirable."
That is the *node-in-cell* predicate (P3), stated as a desideratum, on the **open**
interval — **not** a clamp, and not narrower than (0,1).

### The three angular equations, and the sense in which the quadrature is "fixed"

Printed p. 239 (PDF p. 4), Eqs. (13a)–(13c) — the whole crux of Q1:

| eq. | content | ORPHEUS status |
|---|---|---|
| (13a) | `α_{m+1/2} − α_{m−1/2} = −μ_m w_m` | ✅ same α recursion ORPHEUS uses (Hébert 3.398) |
| (13b) | `(1 − τ_m) α_{m+1/2} + τ_m α_{m−1/2} = (1 − μ_m²)/2` | ⛔ **ORPHEUS does not impose this** |
| (13c) | `μ_m = τ_m μ_{m+1/2} + (1 − τ_m) μ_{m−1/2}` | ✅ = BMC 43, the barycentric τ |

(13b) is not decorative: it is the **first-derivative consistency condition** for
the angular streaming term, obtained by requiring the `∂ψ/∂μ` coefficient in the
Taylor expansion Eq. (7) to reproduce the analytic coefficient `(1−μ²)/r` — see
Eq. (10b) with Eq. (12) substituted.

**Three equations, and the μ-mesh (the cell edges `μ_{m±1/2}`) is given.** So the
unknowns are `α_{m+1/2}`, `μ_m`, and `τ_m` — three unknowns, three equations. The
paper says so explicitly (printed p. 239): the set is to be viewed "as determining
recursively either the three unknowns `α_{m+1/2}`, `μ_m`, and `τ_m` or `α_{m−1/2}`,
`μ_m`, and `τ_m`."

Eliminating `α_{m∓1/2}` and `τ_m` collapses the triple to a **quadratic in the
ordinate** — printed p. 240 (PDF p. 5), Eqs. (14a)/(14b):

`3μ_m² − 2 μ_{m+1/2} μ_m + 2 α_{m−1/2} − 1 = 0`
`3μ_m² − 2 μ_{m−1/2} μ_m + 2 α_{m+1/2} − 1 = 0`

**⭐ THIS is "fixing a quadrature."** Answer to the question as asked:

1. What the scheme takes as **input** is the set of angular **cell widths /
   weights** `w_m` (the abstract: "a set of quadrature weights is specified but
   … the points at which the function is evaluated remain to be determined by the
   analysis," printed p. 238).
2. What it **outputs** is the set of **ordinates** `μ_m` — solved recursively from
   (14a) marching up from `α_{1/2} = 0`, and from (14b) marching down from
   `α_{M+1/2} = 0`. **Yes: choosing the weights forces the ordinates to particular
   values.** They are no longer the Gauss points; they are Eq.-(14) points.
3. The by-product, proved on printed p. 240 as "property 2", is that the resulting
   set satisfies the **diffusion condition automatically**:
   `Σ_m w_m μ_m² = ⅓ Σ_m w_m`, for *any* symmetric set of intervals. The telescoping
   proof needs exactly the seeds `α_{1/2} = α_{M+1/2} = 0`.
4. **And the price, stated by the primary authors themselves** (Conclusions,
   printed p. 248): with a spherical-harmonics treatment of anisotropic scattering,
   "the fourth and higher even angular moments (if present) would not be correctly
   integrated," with the magnitude "likely to be small (<3 % for the fourth moment
   and S₈ Gauss weights) and decreasing for increasing order of quadrature."

Point 4 **is** Lathrop (2000)'s "fixes a quadrature that does not correctly
integrate Legendre polynomials," and here it is in the primary source with a
number attached. The abstract's own framing is the positive half of the same
sentence: "Resulting angular quadrature sets are shown to integrate exactly
polynomials of second order in the direction cosines" — *second* order, and no
more.

### The measured evidence (Table I, printed p. 247 / PDF p. 12)

Six-group homogeneous spherical benchmark (ANL-7416), reference
`k_eff = 0.99597`, 40 uniform spatial intervals. **Percentage error in `k_eff`:**

| S_N | Gauss quad, diamond | Gauss **weights only**, node = interval average, diamond | Gauss weights, μ from Eq. (14), diamond | + space wt only | + **angle wt only** | + space & angle wt | **Equal weights**, μ from Eq. (14), space & angle wt |
|---|---|---|---|---|---|---|---|
| 2 | 3.577 | **13.87** | 3.577 | 3.578 | 2.094 | 2.096 | 2.096 |
| 4 | 0.921 | **4.29** | 0.965 | 0.966 | 0.737 | 0.738 | 0.407 |
| 8 | 0.259 | **1.20** | 0.225 | 0.225 | 0.201 | 0.201 | 0.097 |
| 16 | 0.070 | — | 0.055 | 0.055 | 0.053 | 0.053 | 0.024 |
| 32 | 0.019 | — | 0.014 | 0.014 | 0.014 | 0.014 | 0.006 |

Column 2 is the **δ = 0 ⟹ midpoint-nodes experiment**, run in 1970: keep the
Gaussian weights, but move each ordinate to the centre of its interval so that
τ ≡ ½ is exact. Text, printed p. 246: that set "does not satisfy the diffusion
condition and leads to the poor results shown." **4.6× worse at S₈, 3.9× at S₂.**
This is the primary measured basis for Lathrop (2000)'s δ = 0 argument.

Note also the **role split** the same page states: the *space* weighting (Eq. 11)
"has very little effect on the eigenvalue" but makes "the spurious central flux
dip … disappear"; the *angle* weighting is what "improve[s] the accuracy of the
eigenvalue appreciably." ⚠ In **this** paper the central-flux-dip cure is the
**spatial** weight, not τ. Whether Morel–Montry's "flux dip" is the same object
is not decidable from this source.

**VERDICT:** ANSWERED, and the answer is load-bearing. "Reed and Lathrop's
weighted diamond" = Eqs. (3)/(4) node-weighting **plus** the extra consistency
condition Eq. (13b), whose effect is that the *ordinates become outputs of the
difference scheme* (Eq. 14, a quadratic per ordinate marched from `α_{1/2} = 0`).
It "fixes a quadrature" in the literal sense: the weights are chosen, the points
are then determined, the second moment comes out exact for free, and the 4th and
higher even moments are sacrificed (authors' own estimate <3 % at the 4th moment,
S₈). **It is NOT the same scheme as Morel–Montry's weighted diamond**, which
imposes only (13a)+(13c) and leaves the ordinates alone.

---

## Q2 — Is the ANGULAR truncation error analysed separately from the spatial one?

**Yes, cleanly separated, throughout.** Printed p. 239 (PDF p. 4), Eq. (7) is the
Taylor expansion of the difference equation with the spatial and angular
contributions in two distinct braces. The two are then tracked independently:

- **Spatial coefficient** (of `Δr²`), printed p. 240–241:
  `[(1−a)²A_{i+1} − a²A_i]/V = (5/12)(1/r_{i+1/2}) + O(Δr)` with `a` from Eq. (11),
  giving a second-order error `∝ (Δr²/r)(∂²ψ/∂r²)` — "an improvement over the
  equivalent term `(Δr²/r²)(∂ψ/∂r)` for the standard diamond difference equations."
- **Angular coefficient** (of `w_m²`), printed p. 241 — the one you want.

### The angular error coefficient, and the shape of the amplification

Printed p. 241 (PDF p. 6). The `w_m²` term of Eq. (7) carries the coefficient

`[(1 − τ_m)² α_{m+1/2} − τ_m² α_{m−1/2}] / w_m`

and **property 3** (printed p. 240) is precisely the demand that this **not vary
inversely with `w_m`**. R&L then show that *if* `τ_m = ½ + O(w_m)`,

`[(1 − τ_m)² α_{m+1/2} − τ_m² α_{m−1/2}] / w_m = ¼ [(α_{m+1/2} − α_{m−1/2})/w_m] + O(w_m)/w_m = −μ_m/4 + O(w_m)/w_m`

i.e. the `¼` factors out and the α-difference collapses by Eq. (13a) to `−μ_m w_m`,
leaving a **bounded** coefficient. Their conclusion, verbatim and short: "if
`τ_m = 1/2 + 0(w_m)`, the coefficient of `w_m²` in Eq. (7) is not inversely
proportional to `w_m`, and this term is truly second order."

**⭐ Read the contrapositive.** The *only* thing preventing that coefficient from
being `O(1/w_m)` — i.e. from degrading the angular scheme from second order to
first — is that `τ_m` be within `O(w_m)` of ½. The mechanism is exactly the
`τ²`-vs-`(1−τ)²` **imbalance**: at τ = ½ the two square factors are equal and the
bracket collapses onto the α-*difference* (which is `O(w_m)` by Eq. 13a); away
from ½ the bracket retains an α-*sum*-like remnant that does **not** carry a
factor `w_m`, and the division by `w_m` then bites.

R&L give the deviation in closed form — Eq. (15), printed p. 241:

`τ_m = ½ − (1/(2 w_m)) (μ_{m+1/2} + μ_{m−1/2} − 2 μ_m)`

**i.e. the τ-deviation from ½ is the node's offset from the cell midpoint,
divided by the cell width.** And Eq. (16) is the requirement this places on the
ordinate:

`μ_m = (μ_{m+1/2} + μ_{m−1/2})/2 + O(w_m²)`

**"the ordinate must be the midpoint of its own angular cell to within `O(w²)`."**

⚠ Note the direction of the logic, because it is the opposite of what the phrase
"weighted diamond" suggests: R&L's second-order angular accuracy **requires τ to
be asymptotically ½**, and their whole construction (moving the ordinates) exists
to *make that true* while still satisfying (13b). Their weighted diamond is a
*perturbation of* plain diamond, `τ = ½ + O(w)`, not a departure from it.

### Is there an error-*propagation* factor of the shape `(1−τ)/τ`?

**No.** Searched the full text for amplification/accumulation/stability
vocabulary: the words `amplif`, `accumul`, `stabil`, `oscillat`, `monoton`,
`bound` (in the estimate sense) do not occur anywhere in this paper. There is
**no** von-Neumann-style analysis, no growth factor, no statement about how an
error in `ψ_{m−1/2}` propagates to `ψ_{m+1/2}` through the inverted closure. The
`τ`-dependence R&L expose is entirely a **local truncation-error coefficient**,
not a **recurrence gain**.

The one recursion R&L *do* analyse for propagation is the recursion that
generates the ordinates, not the one that generates the fluxes — see Q4.

**VERDICT:** ANSWERED for the truncation error (yes, separated; the angular error
is second order **iff** `τ_m = ½ + O(w_m)`, Eq. (15)/(16), printed p. 241, and the
governing coefficient is `[(1−τ)²α_{m+1/2} − τ²α_{m−1/2}]/w_m`).
**NOT ADDRESSED for error amplification across the angular sweep** — no `(1−τ)/τ`
gain, no product/recursion over ordinates, no stability analysis of the flux
recurrence anywhere in the paper.

---

## Q3 — Positivity, negative fluxes, flux fix-ups

Grepped the whole text for `negativ|positiv|fixup|fix-up|stabil|oscillat|monoton|
bound|amplif`. Results:

- **The only occurrence of "negative flux" in the paper** is a methodological
  parenthetical in Sec. III, printed p. 246: they modified DTF-IV "(including the
  **removal of negative flux fixup recipes that locally perturb the difference
  schemes**)" so that the comparison would measure the differencing and not the
  fix-up. That is a statement that fix-ups *exist and contaminate a comparison* —
  it is not an analysis, a condition, or a recommendation.
- Every other `negativ`/`positiv` hit is about the **sign of a square root** in
  Eq. (14): "one must choose the negative sign for μ < 0 and the positive sign for
  μ > 0" (printed p. 241). That is root selection in the *ordinate* quadratic, not
  flux positivity.
- The `(0,1)` restriction on `a_{i+1/2}` and `τ_m` (printed p. 238) is justified
  solely by **geometric containment** of the node in its cell ("which seems
  desirable"), with no positivity or stability argument attached.
- ⚠ **One clamp does appear in this paper — but on the SPATIAL weight, and it is
  not theirs.** Footnote 8, printed p. 239: Grant lets his weights depend on the
  sign of `a`, "but this is necessary only to keep `a` between ½ and 1." So the
  `[½, 1]` interval in this literature lineage belongs to **Grant's spatial weight
  `a`**, not to τ. R&L do not adopt it; their own Eq. (11) puts
  `a_{i+1/2} = ½ + (1/6)·Δr/(r_{i+1}+r_i)`, which lies in `[½, ⅔]` by construction
  (⅔ exactly at the innermost cell where `r_i = 0`) — bounded as an *outcome*, not
  imposed. This is worth knowing if the retired ORPHEUS `[½,1]` absorber's
  provenance is ever re-litigated: the interval is real in the literature, and it
  is about `a`, not `τ`.

There is **no** positivity condition, **no** CFL-like or von-Neumann stability
analysis, and **no** statement of when the angular recurrence can drive
`ψ_{m+1/2}` negative — for the spatial scheme, the angular scheme, or either.

**VERDICT: NOT ADDRESSED in this source.** The paper is purely a *local
truncation-error* analysis of a linear difference scheme; positivity is
acknowledged to exist only as the DTF-IV fix-up they deliberately switched off.
No condition here for ORPHEUS to enforce. (Corollary: this paper cannot be cited
as authority that half-angle-flux positivity is un-needed either — it simply
never asks.)

---

## Q4 — The starting direction

**Substantively addressed — but for the ORDINATE recursion, not the flux
recursion.** This is the distinction that decides whether the paper helps you.

### What the paper does say (printed pp. 239–242, PDF pp. 4–7)

- **The seed is `α_{1/2} = 0`, justified by flow direction** (printed p. 239): if
  the M directions "are ordered so that neutrons flow out of, but not into, the
  first μ cell, then `α_{1/2} = 0`. If neutrons flow into, but not out of, the last
  cell, then `α_{M+1/2} = 0`." Both seeds are used; the march for the half-range
  `(−1,0)` runs **forward from `α_{1/2} = 0`** via (14a)+(13a), and the half-range
  `(0,1)` runs **backward from `α_{M+1/2} = 0`** via (14b). Two seeds, two marches,
  meeting at μ = 0.
- **The starting ordinate has a closed form.** Eq. (14a) at m = 1 with
  `μ_{1/2} = −1`, `μ_{3/2} = −1 + w_1`, `α_{1/2} = 0` reduces to
  `3μ_1² − 2(−1 + w_1)μ_1 − 1 = 0`, hence
  `μ_1 = (−1 + w_1)/3 − (1/3)[(−1 + w_1)² + 3]^{1/2}` (negative root), which expands
  to **Eq. (17), printed p. 241**: `μ_1 = −1 + w_1/2 + O(w_1²)`.
  So the starting ordinate is the **midpoint of the first angular cell to
  `O(w₁²)`** — i.e. `τ_1 = ½ + O(w_1)`, same as every other ordinate.
- **⭐ The seed is the base case of an induction that propagates a property
  through the whole level.** Printed pp. 241–242, Eqs. (18)/(19): assuming Eq. (16)
  (`μ_m` = cell midpoint + `O(w²)`) holds at index m, subtracting the two forms of
  Eq. (14) and substituting gives
  `μ_{m+1} = (μ_{m+3/2} + μ_{m+1/2})/2 + O(w_m²) + O(w_{m+1}²)` — Eq. (19) — which is
  Eq. (16) at index m+1, **provided the root sign is chosen correctly**. The paper
  then states that this sign choice "will provide a unique solution to the set of
  Eqs. (13)."

That induction is the closest thing in this paper to your question: it says the
scheme's second-order angular property is **inherited from the m = 1 seed and
carried forward by the recursion**, and that the *wrong branch of the quadratic*
breaks the inheritance. But the quantity being propagated is the **near-midpoint
property of the ordinates**, not an error in ψ.

### What the paper does NOT say

Nothing about how an error in the starting-direction **flux** `ψ_{1/2}` (the
`μ = −1` / `ω = π` seed value) propagates into the bulk ordinates. There is no
seed-flux equation in the paper at all — the `μ = −1` starting-direction balance
equation (the one that has no angular-redistribution term) is never written,
never differenced, and never analysed. Its truncation error is not discussed.

### `[M]` My reproduction of the R&L recursion (not the paper's numbers)

Coded Eqs. (13a)/(14a)/(14b)/(13c) with the paper's two-sided march
(`.venv/bin/python`, this session). Property 2 reproduces exactly
(`Σ w μ² = 0.6666666667` at every order tested, both weight families) and
property 1 (symmetry) holds to `1e-14`, so the implementation is faithful.
Consequences relevant to your amplification question:

| | S8 Gauss wts | S32 Gauss wts | S128 Gauss wts |
|---|---|---|---|
| R&L determined τ range | `[0.4397, 0.5603]` | `[0.4663, 0.5337]` | `[0.4827, 0.5173]` |
| max per-step `\|(1−τ)/τ\|` | 1.274 | 1.144 | 1.070 |
| peak *running* product over the level | 1.574 | 1.705 | 1.772 |
| end-to-end product over the level | **1.000000** | **1.000000** | **1.000000** |

Three things fall out, all of them structural rather than R&L claims:

1. **R&L's own τ never leaves a narrow band around ½**, and the band closes like
   `O(1/N)` exactly as Eq. (15)/(16) requires. Their scheme cannot produce a τ
   small enough to give a large `(1−τ)/τ`.
2. **The end-to-end gain across a symmetric level is pinned to exactly 1**, and
   this is a *theorem, not a coincidence*: property 1 forces the reflection
   identity `τ_{M+1−m} = 1 − τ_m` (verified to `1e-14`), so the factors pair off as
   `[(1−τ)/τ]·[τ/(1−τ)] = 1`. **This holds for any barycentric τ on a symmetric
   partition with symmetric ordinates — including ORPHEUS's.** It is a cheap
   invariant to assert on your own τ table: if your end-to-end product is not 1,
   either the level is not symmetric or the τ is not barycentric on it.
3. The transient therefore has to peak **at the half-level** (μ = 0 crossing) and
   come back down — measured peak 1.5–1.8, growing only logarithmically slowly with
   N. **A measured ~9× transient is an order of magnitude outside what R&L's τ
   produces**, which points the finger at the τ values (or the partition that
   generates them) rather than at the recurrence being intrinsically amplifying.
4. ⚠ `[M]` A *one-sided* march of the α recursion across the whole level (rather
   than R&L's two marches seeded at each end) is ill-conditioned: at S64 with Gauss
   weights it loses the reflection symmetry entirely (deviation 1.18) and the gain
   product blows to `1.07e9`. R&L prescribe the two-sided march for **symmetry**
   (printed p. 240), not stability, but the conditioning consequence is real and is
   worth checking in ORPHEUS's α march if it runs one-sided.

**VERDICT: PARTIALLY ADDRESSED.** The starting *direction* is treated seriously —
`α_{1/2} = 0` seeds the march, Eq. (17) gives the starting ordinate in closed form
(`μ_1 = −1 + w_1/2 + O(w_1²)`), and an explicit induction (Eqs. 18/19) carries the
second-order property from that seed through every subsequent ordinate, with a
root-sign rule as the hinge. **But it is the ordinate/α recursion that is seeded
and propagated, not the flux.** The accuracy of the μ = −1 starting-direction
*flux* equation, and the propagation of *its* error into the bulk ordinates, is
**NOT ADDRESSED**.

---

## Q5 — Cylindrical vs spherical scope

Explicit, and better than Lathrop (2000):

| section | geometry | printed pp. | content |
|---|---|---|---|
| II.A | **1-D spheres** | 238–242 | full derivation, Eqs. (1)–(19) |
| II.B | **3-D cylinders** (r, θ, z, ξ, ω) | 242–243 | full derivation, Eqs. (20)–(28) |
| II.C | **3-D spheres** (r, θ, φ, μ, ω) | 243–246 | full derivation, Eqs. (29)–(45), two angular derivatives |
| III | numerical results | 246–248 | **1-D spherical only** |

**Cartesian is explicitly excluded** (printed p. 238): "Cartesian geometries are
not considered, since no directional derivatives are present in the transport
equation for these geometries."

### The cylinder content, since that is the part you need

Eq. (20), printed p. 242 — the 3-D cylindrical transport equation with the
azimuthal redistribution written as `− (η/r) ∂ψ′/∂ω`, with
`μ = (1−ξ²)^{1/2} cos ω`, `η = (1−ξ²)^{1/2} sin ω`, ω the angle of revolution about
the ξ axis. Change of variable `ψ′(r,θ,z,ξ,ω) → ψ(r,θ,z,ξ,μ,S_η)` (with
`S_η = ±1` tracking the sign of η) uses `∂ψ′/∂ω = −η ∂ψ/∂μ` to give **Eq. (21)**:

`μ ∂ψ/∂r + (η/r) ∂ψ/∂θ + ξ ∂ψ/∂z + (η²/r) ∂ψ/∂μ + σψ = S`

so the cylinder's angular derivative is taken **in the radial cosine μ**, per
ξ-level, exactly as in the sphere — which is why "the case of three-dimensional
cylindrical geometry is only slightly more complicated" (printed p. 242).

The cylinder analogue of the (13) triple is **Eqs. (28a)–(28c)**, printed p. 243:

`α_{m+1/2} − α_{m−1/2} = −μ_m w`
`[(1 − τ) α_{m+1/2} + τ α_{m−1/2}] Δμ = η² w`, with `η² = 1 − μ_m² − ξ_l²`
`μ_m = τ μ_{m+1/2} + (1 − τ) μ_{m−1/2}`

**⭐ Note what (28b) does that the sphere's (13b) hides: it keeps the quadrature
weight `w` and the angular cell width `Δμ` as SEPARATE quantities.** In the sphere
R&L collapse them (`w_m` is defined on printed p. 238 as "quadrature weights,
**normally** `μ_{m+1/2} − μ_{m−1/2}`" — note the hedge), so (13b) has no `Δμ`. The
cylinder derivation never makes that identification: `w` appears in the α
recursion, `Δμ` in the consistency condition, and the telescoping proof of the
diffusion condition (printed p. 243) goes through with them distinct. **That is a
primary-source licence for weight ≠ cell-measure in a cylinder**, which is the
crux my memory flags as "BMC Eq. 52 is NOT a law."

Two more cylinder-specific statements (printed p. 243):

- the seeds are `α_{1/2} = α_{M+1/2} = 0` **at each ξ level**, and "it is not
  necessary that the number of μ-mesh intervals be the same at each ξ level";
- the **ξ levels are free** — "we have produced no restrictions on the choice of
  the ξ levels" beyond symmetry and `Σ ξ² w = ⅓ Σ w`; the per-level μ points then
  satisfy the μ diffusion condition automatically.

The 1-D cylinder is the degenerate case (drop the θ and z terms of Eq. 22); R&L
never write it out, but nothing in §II.B depends on those terms — the angular
system (28) is already θ- and z-free.

**VERDICT: ADDRESSED, and the cylinder IS covered** — §II.B, printed pp. 242–243,
Eqs. (20)–(28), with the angular derivative in the radial cosine per ξ-level.
Contra Lathrop (2000), which has zero cylinder content, this paper carries the
cylinder analogue of every sphere result, including the separated weight-vs-Δμ
form of the consistency condition. **Numerical results are 1-D spherical only** —
no cylinder or 3-D sphere case was ever run.

---

## Q6 — What in this paper CUTS AGAINST what we believe

Six items, ordered by how hard they bite. Two are naming/attribution hazards, two
are genuine cuts, two are supports that come with a caveat.

### C1 ⚠ "Weighted diamond" names two different schemes — ORPHEUS ships the other one

| | Reed & Lathrop 1970 | Morel–Montry / BMC (what ORPHEUS ships) |
|---|---|---|
| imposes (13a) α recursion | yes | yes |
| imposes (13c) barycentric τ | yes | yes |
| imposes **(13b) angular consistency** | **yes** | **no** |
| ordinates `μ_m` | **outputs** (Eq. 14 quadratic) | inputs (from the chosen quadrature) |
| resulting quadrature | exact to degree 2 only; 4th moment off by <3 % at S₈ | whatever the chosen set is |

Not a contradiction of the closure *form* — Eq. (4b) is literally ORPHEUS's
closure. But any ORPHEUS docstring, theory page, or commit message that credits
"Reed & Lathrop's weighted diamond" for the scheme actually implemented is an
L-003-class mis-citation. The safe attribution is Grant (1968) for the node-weight
ansatz Eqs. (3)/(4), Reed & Lathrop (1970) for the ordinate-determining variant,
and Morel–Montry / BMC for the τ-from-a-fixed-quadrature variant.

### C2 ⛔ The genuine cut: R&L's second-order angular accuracy REQUIRES `τ = ½ + O(w)`

Printed p. 241, property 3 + Eqs. (15)/(16). The `w_m²` coefficient in Eq. (7) is
`O(1/w_m)` — i.e. the angular scheme silently degrades from second order to first
— **unless the node sits at its cell midpoint to within `O(w²)`**. R&L's entire
ordinate-moving construction exists to make that true. So this paper is, if
anything, a *pro-midpoint* argument: its weighted diamond is a small perturbation
of plain diamond, not a departure from it.

⚠ Do not over-read it. `τ = ½ + O(w)` is an **asymptotic** statement about the
sequence as `N → ∞`. It does not forbid `τ = 0.42` at a fixed S₈; it forbids
`τ − ½` failing to shrink like `w` as N grows.

⭐ **Which makes it a solve-free gate you can run today.** R&L's own criterion,
restated as a predicate on ORPHEUS's τ table:

> `(τ_m − ½) / w_m` must stay **bounded** as N increases.

`[M]` on R&L's own determined τ this ratio is ≈ ±0.16 at S8-Gauss and stays in the
same band through S128 (τ range narrows `[0.4397,0.5603] → [0.4827,0.5173]` while
`w` narrows proportionally). If ORPHEUS's ratio *grows* with N, R&L's analysis says
the angular differencing is not uniformly second order — and that is a defect in
the **partition/τ rule**, not in the recurrence.

### C3 ⛔ The second genuine cut: their measured accuracy loss lives in the QUADRATURE, not the closure weight and not a recurrence gain

From Table I (printed p. 247), at S₈, decomposing the effects:

- **quadrature second moment**: col 2 → col 3, `1.20 → 0.225` = **5.3×**. Same
  Gaussian weights; the only change is moving the node off the interval average
  onto the Eq. (14) point so the diffusion condition holds.
- **the closure weight τ itself**: col 3 → col 5, `0.225 → 0.201` = **11 %**. Same
  ordinates; the only change is `τ = ½ → τ` from Eq. (13).
- **the spatial weight a**: col 3 → col 4, `0.225 → 0.225` = **nil** on `k_eff`.

So in R&L's configuration the accuracy is dominated by whether the quadrature
integrates `μ²` correctly; the angular closure weight is a ~10 % refinement on top.
Nothing is left over for a recurrence-amplification term to explain.

⭐⭐ **But the reconciliation matters more than the cut.** R&L's configuration
**cannot exhibit your phenomenon**: their τ is pinned near ½ by construction, so
(`[M]`, my reproduction) their per-step gain `|(1−τ)/τ|` never exceeds **1.27** and
their transient never exceeds **1.8**. A ~9× transient is far outside the regime
they analysed. Their silence on amplification is therefore **not evidence of
absence** — it is evidence that in a scheme where τ stays near ½ the question does
not arise. Read the other way round, that is support for your instinct: *if*
ORPHEUS's τ leaves the `½ + O(w)` band, ORPHEUS is operating outside every regime
this literature analysed, and both the truncation error (C2) and the recurrence
gain go bad **together**, driven by the same cause.

### C4 ⚠ Attribution hazard: in THIS paper the flux dip is cured by the SPATIAL weight

Printed p. 246: the space weighting "has very little effect on the eigenvalue.
However, when it is used, the spurious central flux dip … disappears. When the
angular weighting is used, the accuracy of the eigenvalue is improved
appreciably." Two weights, two disjoint payoffs. If any ORPHEUS rationale for the
*angular* weighted diamond rests on "it eliminates the flux dip", this paper
assigns that to `a` (Eq. 11), not to `τ`. Whether Morel–Montry's dip is the same
object (their 1984 TTSP paper is the flux-dip analysis) is **not decidable from
this source** — flagged, not resolved.

### C5 ✅ SUPPORT (with a caveat): τ is determined, not chosen

R&L agree completely that τ is a determined quantity — Eq. (15) writes it in
closed form as the node's offset from the cell midpoint over the cell width. The
caveat is that ORPHEUS's determination uses **one** of R&L's three equations
where R&L use **three**, so the ORPHEUS τ generally leaves a **residual in (13b)**:

> `R_m ≡ (1 − τ_m) α_{m+1/2} + τ_m α_{m−1/2} − (1 − μ_m²)/2`   (sphere, from Eq. 13b)
> `R_m ≡ [(1 − τ_m) α_{m+1/2} + τ_m α_{m−1/2}] Δμ_m − η_m² w_m`   (cylinder, from Eq. 28b)

and by Eq. (7) with Eq. (12) substituted, a non-zero `R_m` leaves the angular
streaming coefficient wrong by **`2 R_m / r_{i+1/2}`**, multiplying `∂ψ/∂μ` — a
*first*-derivative error term that R&L's scheme removes identically. Calibration
for what "acceptable" means: R&L state (printed p. 246) that plain diamond's
residual gives an error `O(Δμ²)·((1−μ²)/2)·∂ψ/∂μ`, i.e. **`R_m` should come out
`O(w²)`**. This is a second solve-free gate, per-ordinate, computable from the τ
table you already have — and unlike the ν-closure gate it is *not* blind on a
symmetric level, because it is a pointwise residual rather than an antisymmetric
sum.

### C6 ✅ SUPPORT: nothing here clamps τ, and the `[½,1]` interval belongs to `a`

R&L restrict `τ ∈ (0,1)` on node-containment grounds only (printed p. 238). The
`[½,1]` interval that appears in this lineage is **Grant's, on the spatial weight
`a`** (footnote 8, printed p. 239), and R&L neither adopt it nor need it — their
Eq. (11) lands `a` in `[½, ⅔]` as an outcome. **Denominator for the standing "no
source clamps τ" claim: this is now 4 of 4 sources that define τ (BMC 2010,
Lathrop 2000, Hébert 2009, Reed & Lathrop 1970) and none clamps it.** The nearest
thing to a clamp in any of them is Grant's `a ∈ [½,1]`, on a different variable.

---

## Notation map (R&L → ORPHEUS), and the conflict to watch

| R&L | meaning | ORPHEUS |
|---|---|---|
| `μ` (sphere, §II.A) | streaming / radial cosine | `mu_x` |
| `τ_m` | angular node weight | `tau` — **same letter, same meaning** |
| `a_{i+1/2}` | spatial node weight | spatial diamond weight |
| `w_m` | angular quadrature weight ("normally `μ_{m+1/2} − μ_{m−1/2}`") | quadrature weight |
| `α_{m±1/2}` | angular redistribution coefficient, `α_{m+1/2} = α_{m−1/2} − μ_m w_m` | `alpha` (= Hébert 3.398, same sign) |
| `α′` | `(A_{i+1} − A_i) α` — the **area-scaled** α | no direct counterpart; do not confuse with `α` |
| `N` | discrete angular flux (Carlson–Lathrop notation) | `psi` |

⛔ **CYLINDER COSINE NAMES ARE A 3-CYCLE AWAY FROM BMC/HÉBERT.** In §II.B R&L use
`μ` = radial, `η` = azimuthal, `ξ` = axial. Bailey–Morel–Chang and Hébert use
`η` = radial, `ξ` = azimuthal, `μ` = axial. Every symbol is occupied by a different
axis in the two conventions:

| axis | ORPHEUS | Reed & Lathrop | BMC / Hébert |
|---|---|---|---|
| radial | `mu_x` | **μ** | η |
| azimuthal | `mu_y` | **η** | ξ |
| axial | `mu_z` | **ξ** | μ |

So R&L's cylinder Eq. (28b) `η² w` is ORPHEUS's *azimuthal* cosine squared, and
R&L's `μ_m` in the same equation is ORPHEUS's *radial* cosine. Transcribing (28)
symbol-for-symbol against a BMC-conditioned reading gives the wrong equation.

---

## Provenance ledger

| claim | status |
|---|---|
| Reed & Lathrop, NSE **41**(2) 237–248, 1970-08, DOI `10.13182/NSE70-A20710` | ✅ **CrossRef-verified** this session — matches the brief exactly |
| every equation quoted (Eqs. 2, 3, 4, 7, 9–17, 19, 20–25, 28) | ✅ read on the rendered scan, PDF pp. 3–7, not sidecar-only |
| Table I values | ✅ from `.mocr.json` `pages[11].tables[0]`; header structure re-read on PDF p. 12 |
| `printed = PDF + 235` | ✅ confirmed on rendered running heads (238…248) |
| the τ / gain / `Σwμ²` numbers marked `[M]` | ⚠ **my reproduction of R&L's recursion, not the paper's numbers** — labelled as such at each occurrence; property 1 and property 2 reproduce to `1e-14`, which is the evidence the implementation is faithful |
| Grant, "Numerical analysis of discrete ordinate methods," *J. Comput. Phys.* **2**(4) 381–402 (1968), DOI `10.1016/0021-9991(68)90044-2` | ✅ **CrossRef-verified** — R&L footnotes 6 and 8 cite it correctly |
| Morel–Montry 1984 "flux dip" = the dip R&L cure with the spatial weight | ⛔ **NOT decidable from this source** — flagged in C4, not resolved |
| Carlson & Lathrop 1968 (their ref. 5) and Lathrop & Carlson JCP **1**, 173 (1966) (their ref. 7) | not chased; cited here only as R&L's own attributions |

### ⭐ Acquisition lead — Grant (1968) is the primary this paper builds on, and it is NOT local

`scratch/literature/` has no Grant. Yet R&L's Eqs. (3)/(4) — the entire
weighted-diamond node ansatz, i.e. the closure ORPHEUS ships — are attributed in
the text to "a suggestion due to Grant" (footnote 6, printed p. 238), and footnote
8 (printed p. 239) says Grant "does not determine angular weights and limits his
degrees of freedom by assuming an explicit form for α," while needing a
sign-dependent weight "only to keep `a` between ½ and 1."

So Grant 1968 is (a) the origin of the weighted-diamond ansatz, (b) the origin of
the `[½,1]` weight interval that ORPHEUS recently retired, and (c) an *alternative*
determination of α — three things directly on the live question, in a paper nobody
in this corpus has read. **This is the L-012 pattern exactly**: the local text names
its own primary, and that primary is unread.

⚠ OA status checked (L-009): `is_oa = False`, `oa_status = closed`, no repository
copy. Elsevier JCP 1968 — a real acquisition decision, not a Cloudflare 403. Your
call whether to add it.

## The one-line answer

**Yes — but not the way you expect: Reed & Lathrop never analyse error propagation
through the angular recurrence, and instead supply the two solve-free criteria
that decide whether your τ can ever get into the amplifying regime at all.** Their
scheme pins `τ = ½ + O(w)` as a *precondition* for second-order angular accuracy
(Eq. 15/16, printed p. 241) and, `[M]` reproduced here, their own determined τ
never leaves `[0.44, 0.56]` at S₈, capping the per-step gain at 1.27 and the
transient at 1.8 with the end-to-end product pinned to exactly 1 by the reflection
identity `τ_{M+1−m} = 1 − τ_m`. A measured ~9× transient is therefore an order of
magnitude outside anything this literature contemplates, which means the τ values
themselves — not the recurrence's existence — are the anomaly, and the two gates
to run before touching the seed are **(i)** `(τ_m − ½)/w_m` bounded under N-refinement
and **(ii)** the Eq.-(13b)/(28b) residual `R_m` scaling as `O(w²)`. On the seed
itself the paper is silent about flux error, but it does give `μ_1 = −1 + w_1/2 +
O(w_1²)` (Eq. 17) and an explicit induction (Eqs. 18/19) proving the *ordinate*
property propagates from the `α_{1/2} = 0` seed — with a root-sign rule as the
hinge, and, `[M]`, with the α march ill-conditioned if it is run one-sided across
the whole level instead of R&L's two marches seeded at each end.
