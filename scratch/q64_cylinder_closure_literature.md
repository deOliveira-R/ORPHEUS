# Q6.4 — what the literature prescribes for the CYLINDER's AZIMUTHAL angular-redistribution closure

Extraction date 2026-08-11. Sources, all **LOCAL** in `scratch/literature/`:

| tag | file | printed-page map |
|---|---|---|
| **H** | `Hebert(2009)Chapter3.pdf` (*Applied Reactor Physics*, ch. 3) | printed = PDF + 66 `[VERIFIED p.72→138, p.73→139]` |
| **BMC** | `Bailey-Morel-Chang(2010)Asymptotic Diffusion-Limit Accuracy of Sn Angular Differencing Schemes.pdf` (NSE 165(2):149-169) | printed = PDF + 147 |
| **L** | `Lathropp(2000)A Comparison of Angular Difference Schemes for One-Dimensional Spherical Geometry SN Equations.pdf` (NSE 134:239-264) | printed = PDF + 237 |

Nothing had to be acquired; no pivot to a secondary source was needed.
**Zotero not consulted** (prior sessions this week found the server down; no
annotations checked — see the note in §7).

## ⛔ Two premises in the brief are wrong, and they matter

1. **Hébert's cylinder is §3.9.3, NOT §3.9.4.** §3.9.4 (PDF pp. 76-80, printed
   142-146) is the **1-D SPHERE**. The cylinder is **§3.9.3, PDF pp. 71-75 =
   printed 137-141**, Eqs. **(3.390)-(3.417)**. Every equation in the
   3.418-3.439 range the brief quotes is spherical.
2. **(3.437)/(3.439) are NOT "a weighted one".** They are the *same* τ ≡ ½
   relation as (3.431), algebraically re-arranged into an extrapolation
   `φ_{n+1/2,i} = 2φ_{n,i} − φ_{n−1/2,i}`, and they are the SPHERE's polar-index
   pair. The cylinder's azimuthal counterparts are **(3.412)** (inward spatial
   sweep) and **(3.414)** (outward). No weighted-diamond appears anywhere in
   Hébert ch. 3, in either geometry.

---

## 1. Hébert §3.9.3 — what closure does he prescribe for the CYLINDER's half-angle azimuthal flux?

**Source**: H, §3.9.3, **PDF p. 73 = printed p. 139**, Eq. **(3.406)**
`[VERIFIED-ON-SCAN]`.

The single equation carries **both** closures, spatial and angular, as a chain
of two equalities:

```
φ_{p,q,i} = ½(φ_{p,q,i−1/2} + φ_{p,q,i+1/2}) = ½(φ_{p,q−1/2,i} + φ_{p,q+1/2,i})
             └──── spatial diamond (index i) ───┘   └── ANGULAR diamond (azimuthal index q) ──┘
```

He labels it *"the diamond differencing scheme"* and presents it in the running
text as **the** discretization — not as one option in a comparison, and not as an
ablation. There is no alternative angular closure offered anywhere in §3.9.3.

The recursion actually used in the sweep is the rearranged form, **PDF p. 74 =
printed p. 140, Eqs. (3.412) and (3.414)**:
`φ_{p,q+1/2,i} = 2φ_{p,q,i} − φ_{p,q−1/2,i}` — *identical* in both, i.e. the
azimuthal half-flux extrapolation is the same for µ<0 and µ>0; only the spatial
half of the pair flips (`φ_{p,q,i−1/2}` vs `φ_{p,q,i+1/2}`).

**Which index is which** (this is the brief's actual question):

| Hébert eq. | geometry | the ANGULAR index it closes | variable of that index |
|---|---|---|---|
| **(3.406)** | **cylinder** | `q` = **azimuthal**, base points `ω_{p,q}` on `0 ≤ ω ≤ π` within polar level `ξ_p` | **the angle ω** |
| (3.412)/(3.414) | cylinder | same `q` (extrapolated form) | the angle ω |
| (3.431) | sphere | `n` = polar, base points `µ_n` | the cosine µ |
| (3.437)/(3.439) | sphere | same `n` (extrapolated form) | the cosine µ |

**VERDICT: Hébert prescribes the PLAIN ANGULAR DIAMOND, τ ≡ ½, for the
cylinder's azimuthal half-flux — Eq. (3.406) (printed p. 139), sweep form
(3.412)/(3.414) (printed p. 140). It is his default and only scheme, not an
ablation. τ ≡ ½ is a pure arithmetic mean of the two half-angle fluxes, so it is
VARIABLE-AGNOSTIC: the midpoint of `[½,½]` weights is the same number whether
one thinks of the cell as an interval in ω or in η.**

---

## 2. Does Hébert define a weight τ at all for the azimuthal index?

**No.** The letters τ, δ and the phrase "weighted diamond" do not occur in
§3.9.3 or §3.9.4. Grep over the whole chapter finds no angular weight parameter
in either geometry.

What he *does* define for the azimuthal index is the **α / half-cosine
recursion**, and it is worth being exact about what it is and is not, because it
is the object most easily mistaken for a τ.

**Source**: H, §3.9.3, **PDF p. 72 = printed p. 138**, Eqs. **(3.393)**,
**(3.396)**, **(3.397)**, **(3.398)**, **(3.399)** `[VERIFIED-ON-SCAN]`.
Attributed in the text to **Alcouffe and O'Dell**, his ref. **[36]** (that
reference is NOT in the local folder).

- **(3.393)** — the redistribution integral over the angular cell is evaluated
  as `W_p [η_{p,q+1/2} φ(ρ,ξ_p,ω_{p,q+1/2}) − η_{p,q−1/2} φ(ρ,ξ_p,ω_{p,q−1/2})]`.
  Note what the differencing variable is: the ω-derivative `∂/∂ω[ηφ]`
  integrates **exactly** in **ω** by the fundamental theorem of calculus, so the
  half-fluxes live at **ω-cell edges `ω_{p,q±1/2}`**, and each is multiplied by a
  coefficient labelled `η_{p,q±1/2}`.
- **(3.396)** — the seed: `η_{p,2F(p)+1/2} = 0`, justified geometrically
  ("the first value η lies on the ξ-µ plane and is therefore zero"). This is the
  chapter's *only* geometric anchor for a half-cosine.
- **(3.397)/(3.398)** — every *other* `η_{p,q±1/2}` is obtained **not**
  geometrically but from **preservation of a constant flux in an infinite
  non-absorbing source-free medium**:
  `η_{p,q+1/2} = η_{p,q−1/2} + (W_{p,q}/W_p) · µ_{p,q}`.
- **(3.399)** — with `α_{p,q±1/2} ≡ W_p η_{p,q±1/2}`, this is
  `α_{p,q+1/2} = α_{p,q−1/2} + W_{p,q} µ_{p,q}`.

**Shape, in the brief's requested terms.** The recursion is a **running sum of
(weight × RADIAL cosine µ)** that lands in the **AZIMUTHAL cosine η**. There is
no numerator/denominator difference quotient at all — it is not an interpolation
weight, it is a first-moment conservation condition. Written as a difference it
reads `η_{p,q+1/2} − η_{p,q−1/2} = (W_{p,q}/W_p) µ_{p,q}`: a difference in the
**azimuthal direction cosine** on the left, and on the right the **radial
direction cosine** — *no angle appears anywhere in it*, and no ω-difference
appears anywhere in it.

⚠ **Consequence, and it is the crux of the ω-vs-η question for Hébert.** The
symbols `ω_{p,q±1/2}` appear in (3.393) **only as labels for where the
half-angle flux lives**. No equation in the chapter ever assigns them a
numerical value, and `η_{p,q±1/2}` is *not* `sinθ_p sin ω_{p,q±1/2}` for any ω —
it is whatever (3.398) makes it. So Hébert's azimuthal cell has a
**recursion-defined edge in the cosine and an unevaluated label in the angle**,
and the diamond closure (3.406) never needs either: `½ + ½` requires no cell
measure.

**VERDICT: NOT ADDRESSED — Hébert defines NO τ for the azimuthal index (nor for
the polar one). His azimuthal machinery is the α/η half-cosine recursion
(3.396)-(3.399), whose differences are in the DIRECTION COSINE (azimuthal η
incremented by weight × radial µ), never in the angle ω. Because his closure is
τ ≡ ½, the question "in which variable is the weight taken?" does not arise in
his scheme at all.**

---

## 3. BMC 2010 — is the CYLINDER's τ defined in the same variable as the sphere's?

**Yes, identically, and they say so in the sentence that introduces it.**

### 3a. The sphere (§IV) — Eqs. (42)/(43)
**Source**: BMC, **PDF p. 8 = printed p. 155** `[VERIFIED-ON-SCAN]`.

- **(42)** `τ_m = (µ_m − µ_{m−1/2}) / (µ_{m+1/2} − µ_{m−1/2})`
- **(43)** `µ_m = τ_m µ_{m+1/2} + (1 − τ_m) µ_{m−1/2}`, introduced as *"uniquely
  defined by the weighted diamond relationship between the cell-center and
  cell-edge **cosines**"*, and (43) *"implies that Eq. (15) will exactly relate
  the cell-edge and cell-center fluxes when the angular flux assumes the linear
  form"* — i.e. τ is **the barycentric coordinate of the node cosine between its
  two edge cosines**, and it is exact for ψ affine in µ. `[CONFIRMED]`

### 3b. The cylinder (§V, R-Z) — Eq. (74)
**Source**: BMC, **PDF p. 13 = printed p. 160**, Eq. **(74)** `[VERIFIED-ON-SCAN]`.

```
τ_{m,n} = (µ_{m,n} − µ_{m−1/2,n}) / (µ_{m+1/2,n} − µ_{m−1/2,n})
```

introduced verbatim as *"the Morel and Montry weighting factors in Eq. (53)
analogous to those defined by Eq. (42)"*. In BMC's R-Z convention (Fig. 2,
printed p. 156) **µ = sin θ cos ω is the RADIAL cosine** and `n` indexes the
ξ-levels (ξ = cos θ, axial). So the cylinder's τ is the barycentric coordinate of
the ordinate's **radial cosine** inside a **radial-cosine cell**, computed
independently on each ξ-level. **The azimuthal angle ω does not appear in Eq.
(74), nor in Eq. (53), nor in Eq. (75).**

### 3c. Eq. (52) — where the cylinder's cell EDGES come from
**Source**: BMC, **PDF p. 10 = printed p. 157** `[VERIFIED-ON-SCAN]`.

Prefaced by *"The cell-edge cosines … are obtained as in the spherical-geometry
case except that there is a separate recursion and starting cosine for each
ξ-level"*:

- recursion `µ_{m+1/2,n} = µ_{m−1/2,n} + w̄_{m,n}` (cumulative **weight**),
- seed `µ_{1/2,n} = −√(1 − ξ_n²)` (= sinθ_n·cos π, the inward starting direction),
- terminus `µ_{M+1/2,n} = +√(1 − ξ_n²)`, which *"will always"* be reached when
  computed recursively,
- **weights renormalised per level to `Σ_m w̄_{m,n} = 2√(1 − ξ_n²)`**.

⛔ **PUBLISHED TYPO, verified on the rendered page (L-010 class).** The printed
recursion reads `µ_{m+1/2,n} = µ_{m+1/2,n} + w̄_{m,n}` — self-referential, hence
vacuous. The RHS subscript must be `m−1/2`. Proven three ways: the sphere's
analogous **Eq. (12)** is printed correctly, the stated terminus cannot otherwise
be reached, and the seed would be pointless. **Eq. (50)** (the α recursion) has
the *same* typo on the *same* page — printed `α_{m+1/2,n} = α_{m+1/2,n} −
µ_{m,n} w_{m,n}`, must be `α_{m−1/2,n}`.

### 3d. Eq. (53) — the closure itself, and its stated admissible range
**Source**: BMC, printed p. 157 `[VERIFIED-ON-SCAN]`.

`ψ_{m,n} = τ_{m,n} ψ_{m+1/2,n} + (1 − τ_{m,n}) ψ_{m−1/2,n}` — *"a general
weighted diamond relationship between the cell-edge and cell-center cosines on
each ξ-level"*, with *"the weighting factors can take on **any value between zero
and one**"*; `τ = 1` → step, `τ = ½` → diamond. Eq. **(54)** is the inverted
recursion `ψ_{m+1/2,n} = ψ_{m,n}/τ − ((1−τ)/τ)ψ_{m−1/2,n}`.

### 3e. Eq. (75) — the cylindrical β
**Source**: BMC, printed p. 160 `[VERIFIED-ON-SCAN]`.

`β_cylindrical = Σ_n Σ_m µ_m [α_{m+1/2,n} µ_{m+1/2,n} − α_{m−1/2,n} µ_{m−1/2,n}]`
(the sphere's is the same sum over one level with a ½ prefactor). **Again every
symbol is a radial cosine µ; no η, no ω.** Table II (printed p. 160): MM gives
`−4.12E−16 … −1.55E−15` at orders 2/4/8/12/16, diamond gives
`1.42E+00 / 1.15E−01 / 2.09E−02 / 1.12E−02 / 8.82E−03`.

⚠ **Do not read the `µ_{m±1/2,n}` in (75) as the Eq.-(52) cumulative-weight
edges.** `[M]` 2026-08-11 (prior session, `scratch/q64_tau_edge_convention_literature.md`):
under that reading β is identically zero for *every* scheme, contradicting BMC's
own Tables I/II. Reproducing Table I to every printed digit required the
**scheme-IMPLIED** edge march — Eq. (43) solved forward from the τ under test.
The two objects share a symbol.

**VERDICT: The cylinder's τ IS defined in the same variable as the sphere's —
a barycentric coordinate in the RADIAL DIRECTION COSINE µ (BMC Eq. 74 ≡ Eq. 42,
per ξ-level), with cell edges from the cumulative-WEIGHT recursion Eq. (52) and
per-level weight renormalisation to `2√(1−ξ_n²)`. BMC never poses anything in
the azimuthal angle ω.**

---

## 4. Lathrop 2000 — is τ ≡ ½ a DEFICIENCY or the REFERENCE case?

⚠ **First, a scope fact that bears on how much weight this source can carry
here: Lathrop 2000 is SPHERE-ONLY.** The string "cylind" occurs **zero** times
in the whole paper (title: *"…for One-Dimensional Spherical Geometry S_N
Equations"*). Any cylinder use of it is by analogy, mediated by BMC §V.

His δ is τ affinely: **Eq. (23)**, printed **p. 245** (PDF p. 8)
`[VERIFIED-ON-SCAN]`:
`µ_m = ((1+δ)/2)µ_{m+1/2} + ((1−δ)/2)µ_{m−1/2} = µ̄ + δΔµ/2`, i.e.
**τ = (1+δ)/2**, so **δ = 0 ⟺ τ = ½ ⟺ µ_m = µ̄ (the node is the cell midpoint)**.
Eq. (23) ≡ BMC Eq. (43); Eq. (24) ≡ BMC Eq. (15)/(53).

**The argument runs against τ ≡ ½, in both places, for a GIVEN quadrature.**

1. **Printed p. 245** (just above Eq. 23) `[VERIFIED-ON-SCAN]` — the diamond is
   O(Δµ²) *only* when the node is the midpoint; for an arbitrary quadrature it is
   not, *"both because of the alpha coefficients and because the simple diamond
   scheme corresponds to `ψ_m = ψ(r, µ̄)`, not `ψ_m = ψ(r, µ_m)`."* He then
   introduces the weighted diamond Eq. (24) as *"an appropriate"* remedy, and
   quantifies the diamond's error as **O(δΔµ)** — small for Gauss-Legendre
   *"except near the ends of the µ range"*.
2. **Printed pp. 249-250** (the sentence straddles the page break)
   `[VERIFIED-ON-SCAN]` — after re-deriving Morel & Montry's analysis: *"the
   standard angular diamond scheme, δ = 0, implies a cell-centered µ_m = µ̄, which
   … would not have c = ⅔ and hence would not give the correct diffusion limit"*
   (`c ≡ Σ µ_m² Δµ_m`, **Eq. 54**; the diffusion condition `c = ⅔` is **Eq. 29**).

⭐ **The apparent near-opposition in the brief dissolves once "δ = 0 as a
property of the NODES" is separated from "τ ≡ ½ as a CLOSURE imposed on given
nodes".** The p. 245 sentence *"only with µ_m = µ̄ (δ = 0) is the truncation order
O(Δµ²)"* sits inside §III.C, which is about **Reed & Lathrop choosing the
cosines**; it says *pick a quadrature whose nodes are the midpoints*. Printed
p. 250 then closes that door quantitatively: **Eq. (62)** gives the equal-interval
midpoint nodes `µ_m = −1 + (2m−1)Δµ/2` and **Eq. (63)** their
`c = (2/3)(1 − 1/N²) = ⅔ − Δµ²/6` — so the δ = 0 quadrature *provably* misses the
diffusion condition, by `O(Δµ²)`. `[M]` Eq. (63) reproduces the previously
measured midpoint-c ladder exactly: N = 2/4/8/16 → 0.5 / 0.625 / 0.65625 /
0.6640625.

Lathrop's own numerics take the midpoint *"except for the weighted diamond
scheme, where the direction µ_m is calculated as part of the approximation"*
(printed p. 250), and his §VII conclusion ranks the **weighted** diamond above
the plain diamond: *"the weighted diamond difference approximation has smaller
maximum and average relative flux errors than the diamond"* (abstract, printed
p. 240; conclusion §VII, printed p. 262 — *"the weighted diamond difference
approximation has the smallest initial error"*).

**VERDICT: Lathrop treats τ ≡ ½ (δ = 0) as a DEFICIENCY, not as the reference or
exact case — twice, on two independent grounds (a truncation argument on printed
p. 245: the diamond returns ψ at µ̄ rather than at µ_m, error O(δΔµ); and a
diffusion-limit argument on printed pp. 249-250 with the closed form Eq. 63:
δ ≡ 0 forces midpoint nodes, whose `c = ⅔(1−1/N²) ≠ ⅔`). His weighted δ from
Eq. (23) is the prescribed remedy. ⚠ The verdict is SPHERICAL — he never treats
the cylinder.**

---

## 5. Does any source prescribe a LIMITER / CLAMP on τ (e.g. τ ∈ [½,1])?

**No. Not one of the three, in either geometry. And BMC's own worked example
violates a τ ≥ ½ rule.**

### 5a. What the sources DO state as the admissible range — always `[0,1]`

| source | statement | where |
|---|---|---|
| BMC (sphere) | *"each weighting factor τ_m can take on any value between zero and one"* | printed p. 152, at Eq. (15) |
| BMC (cylinder) | *"the weighting factors can take on any value between zero and one"* | printed p. 157, at Eq. (53) `[VERIFIED-ON-SCAN]` |
| Lathrop | `|δ| ≤ 1` *"as long as µ_m is in Δµ"* ⟺ `τ ∈ [0,1]` | printed p. 245, at Eq. (23) `[VERIFIED-ON-SCAN]` |

Note Lathrop's phrasing: `τ ∈ [0,1]` is a **consequence** of the ordinate lying
inside its own cell, not an imposition. It is a well-posedness *predicate*, and
the thing to test is the premise (node ∈ cell), not the conclusion.

### 5b. Candidate justification (a) — positivity / non-negative flux: NOT ADDRESSED

The only positivity analysis in the corpus is **SPATIAL**: Hébert §3.9.2,
printed p. 137, Eqs. **(3.387)-(3.389)** — the slab diamond's mesh-averaged flux
goes negative unless `Δx_i Σ_i < 2|µ_n|`. There is **no angular analogue**
anywhere in ch. 3, in any geometry. Lathrop mentions negative fluxes exactly
once, explicitly as a *"foible … of a particular **spatial** difference
approximation"* (printed p. 262), and his one strictly-positive citation
(ref. 2, Walters & Wareing 1996) is a spatial characteristic scheme. BMC: silent.

### 5c. Candidate justification (b) — stability of the angular recurrence: NOT ADDRESSED

BMC **do** write the amplifying recurrence — Eq. (54), printed p. 157:
`ψ_{m+1/2,n} = ψ_{m,n}/τ − ((1−τ)/τ)ψ_{m−1/2,n}` — and they say nothing about
`|(1−τ)/τ| ≤ 1`, about error amplification, or about τ ≥ ½. Neither does
Lathrop (his Eq. 31 recasts γ_m with the factor `(1−δ)/(1+δ)` = `(1−τ)/τ`, again
with no stability remark). Neither does Hébert. **The `τ ≥ ½ ⟺ |(1−τ)/τ| ≤ 1`
argument is not in any of the three sources.**

### 5d. The refutation of a universal τ ≥ ½ — BMC Eq. (47) CONFIRMED

**Source**: BMC, printed **p. 155** (PDF p. 8), Eq. **(47)** `[VERIFIED-ON-SCAN]`,
in the S₂ **spherical** worked example. The printed pair is `τ₁ = µ₁ + 1` and
`τ₂ = µ₂`, and **Fig. 1 on the same page labels the S₂ set explicitly**:
starting direction `µ_{1/2} = −1`, `µ₁ = −(1/3)^{1/2}`, interpolated
`µ_{3/2}`, `µ₂ = (1/3)^{1/2}`, ending `µ_{5/2} = 1`. Substituting:

```
τ₁ = 1 − 1/√3 = 0.42265…  < ½        τ₂ = 1/√3 = 0.57735…  > ½
```

`[M]` reproduced from BMC's own construction (Eq. 12 cumulative-weight edges on
Gauss-Legendre + Eq. 42): S₂ → `[0.42265, 0.57735]`, matching Eq. (47) to all
printed digits. **The brief's number and its context are CONFIRMED.** The
sub-½ half is not an S₂ curiosity — `[M]` at every order exactly **half** the set
is below ½ (S₄ 2/4, S₈ 4/8 → `0.3923, 0.4591, 0.4809, 0.4942`, S₁₆ 8/16), and
Eq. (46) `µ_{3/2} = 0` is *derived from β = 0*, i.e. the sub-½ value is what the
diffusion-limit requirement produces.

BMC also close the door on the inverse move (imposing τ ≡ ½ and letting the
quadrature follow), printed p. 155: a simple average *"would force the quadrature
set to be the midpoint rule, µ₁ = −½, µ₂ = ½. However, the midpoint rule is not a
good choice of quadrature sets in general."* — the same argument Lathrop makes
with Eq. (63).

**VERDICT: NOT ADDRESSED — no source prescribes any clamp or limiter on τ; all
three state the admissible range as `[0,1]` (Lathrop derives it, as a consequence
of node ∈ cell). Neither the positivity nor the stability justification appears
anywhere for the ANGULAR recurrence (positivity is analysed only for the SPATIAL
diamond, Hébert 3.387-3.389). A universal `τ ≥ ½` is REFUTED by BMC's own Eq.
(47): τ₁ = 1 − 1/√3 ≈ 0.4226 for S₂ Gauss-Legendre, CONFIRMED on the printed
page and reproduced numerically. A `[½,1]` floor would corrupt half of every
Gauss-Legendre τ set.**

---

## 6. ⭐ The ω-vs-cosine choice — is it discussed anywhere?

### 6a. Explicit discussion: NONE, in any of the three

- **Lathrop 2000**: no cylinder at all (zero occurrences of "cylind").
- **BMC 2010**: ω appears in exactly two places — the continuous R-Z equation
  **Eq. (48)** (printed p. 156), whose redistribution term is
  `−(1/r)·∂/∂ω[η ψ]`, i.e. **a derivative in the ANGLE ω**; and the coordinate
  definition on the same page. From Eq. (49) onward ω is gone: the discrete
  redistribution is `(α_{m+1/2,n}ψ_{m+1/2,n} − α_{m−1/2,n}ψ_{m−1/2,n})/(r w_{m,n})`,
  and the closure weight Eq. (74) is barycentric in the **radial cosine µ**. They
  never remark on the switch, never define an ω-cell edge, and never mention the
  azimuthal cell measure.
  ⛔ **Published typo on printed p. 156** `[VERIFIED-ON-SCAN]`: the coordinate
  line reads *"ξ = cos θ, µ = sin θ cos ω, η = sin θ cos ω"* — the second must be
  `sin θ sin ω`, else η ≡ µ and `Ω = µ e_r + ξ e_z + η e_γ` degenerates. Fig. 2
  on the same page shows ω measured from `e_r`, and Eq. (52)'s seed
  `µ_{1/2,n} = −√(1−ξ_n²)` (= sinθ·cos π) confirms `µ = sinθ cos ω`. L-010 class.
- **Hébert 2009**: writes the redistribution as `∂/∂ω[ηφ]` and integrates it
  **exactly in ω** — Eq. (3.393), printed p. 138 — so his half-fluxes are
  explicitly at ω-edges `ω_{p,q±1/2}`. But he **never assigns those edges a
  value**, never states the boundaries of the *"angular sub-domain surrounding
  direction Ω_n"* (printed p. 137), and never compares the two variables. His
  closure (3.406) is τ ≡ ½, which needs no cell measure in either variable.

### 6b. ⭐⭐ The highest-value finding: Hébert's own cylindrical quadrature is EQUISPACED IN ω WITH UNIFORM WEIGHTS, so his τ ≡ ½ *is* the barycentric-in-ω weight

**Source**: H, §3.9.1, printed **p. 133** (recommendation) and printed **p. 135**
(the rules), Eqs. **(3.369)-(3.378)** `[VERIFIED-ON-SCAN p. 69 = printed 135]`.

- Printed p. 133: for 1-D cylindrical geometry *"it is recommended to use the
  combination of a Gauss-Legendre and Gauss-Chebyshev quadrature"* with ξ along
  the cylinder axis and µ along the radius.
- **(3.370)** (triangular set, `M ≡ N−2p+2` points on level p) and **(3.375)**
  (product set, `M ≡ N`) — `[M]` both simplify **exactly** to

  ```
  ω_q = π (q − ½) / M          [verified to ≤4.4e-16 for N = 4, 8, 12, 16, 32]
  ```

  i.e. **the MIDPOINTS of M equal ω-cells of width Δω = π/M**.
- **(3.371)/(3.376)** — `W_{p,q} = πW_p/M`, **uniform in q**: the weight is the
  **ω-measure** of that cell (up to the level factor), not its cosine measure.

⟹ **On Hébert's own recommended cylindrical rule the node is the ω-barycentre of
its own cell by construction.** So his `τ ≡ ½` (3.406) is *not* an arbitrary
"plain diamond" on that rule — it is exactly BMC Eq. (43) / Lathrop Eq. (23)
**posed in ω**, i.e. δ = 0 *because the nodes really are the midpoints in the
march variable*. (This is an INFERENCE from his quadrature + closure; he does not
say it.)

### 6c. ⭐ The measure claim in the brief is CORRECT, and it is what breaks the transplant

`[M]` On Hébert's equispaced-ω rule, the **radial-cosine measure of an equal-ω
cell** is `|Δµ| = 2 sinθ · sin(ω_q) · sin(Δω/2)` — **∝ sin ω_q, not constant** —
while the quadrature weight is constant in q. Ratio of the two, per level:

| M | µ-measure(cell) / uniform-µ-width | spread across the level |
|---|---|---|
| 8 | 0.3045 … 1.5307 | **5.0×** |
| 16 | 0.1537 … 1.5607 | **10.2×** |
| 32 | 0.0770 … 1.5683 | **20.4×** |

(the spread ≈ `1/sin(π/2M)` ≈ `2M/π`, so it **worsens with refinement**).

⟹ **BMC Eq. (52) is not a law; it is the statement that in THEIR quadrature the
weight equals the cell's µ-measure.** Consequence, `[M]` on Hébert's rule with
Eq. (52) edges + Eq. (74):

| M | τ from BMC (74) on (52) edges | inside `[0,1]`? |
|---|---|---|
| 8 | `−0.3259 … +1.3259` | **NO** |
| 16 | `−1.1841 … +2.1841` | **NO** |
| 32 | `−2.8552 … +3.8552` | **NO** |

— it violates Lathrop's `|δ| ≤ 1` premise (node ∈ cell) and diverges with
refinement. For contrast, on the **sphere** with Gauss-Legendre the identical
construction stays put: `[M]` τ ∈ `[0.4226, 0.5774]` (S₂) … `[0.3904, 0.6096]`
(S₁₆), terminus `+1.000000000000000` exactly.

Two well-posed readings therefore remain on an equispaced-ω rule, and `[M]` they
are numerically close:

| reading | τ | bounded? |
|---|---|---|
| barycentric in **ω** (arc cell, node = midpoint) | **exactly ½** | trivially |
| barycentric in **µ** on the **geometric arc edges** `sinθ cos(ω_q ∓ Δω/2)` | `½ + ½ cot(ω_q) tan(Δω/4)` (verified to 6e-16) | `[¼, ¾]` at every M, → ½ |

**VERDICT: NOT ADDRESSED explicitly by ANY of the three sources — none of them
remarks on the azimuthal cell's non-uniform cosine measure, and none discusses
which variable the azimuthal closure should be posed in. But the strongest
available IMPLICIT evidence points at ω: Hébert's recommended cylindrical
quadrature (3.370)/(3.375) places every ordinate exactly at the MIDPOINT IN ω of
an equal-ω cell and gives it the cell's ω-measure as its weight
(3.371)/(3.376) — `[M]` verified to 4e-16 — so his τ ≡ ½ (3.406) IS the
barycentric/Morel-Montry weight taken in ω. Conversely the brief's measure claim
is `[M]` CONFIRMED (µ-measure ∝ sin ω, spreading 5×/10×/20× at M = 8/16/32),
and it is exactly why BMC's cosine-partition construction (Eq. 52) cannot be
transplanted onto an equispaced-ω rule: τ leaves `[0,1]` and diverges with
refinement.**

---

## 7. Provenance and gaps

- Everything above is from the **local** folder; nothing was acquired online, and
  **no pivot to a secondary source was made**.
- **Zotero: not consulted — no user annotations checked.** (Prior sessions this
  week recorded the server down; the delegation rule says record and proceed.)
- **NOT LOCAL, and each would sharpen a verdict:**
  1. **Morel & Montry (1984), TTSP 13(5):615** — the primary source for the
     weighted τ. BMC and Lathrop both only *cite* it. It is the one place a
     clamp/positivity argument could still be hiding, and the only place the
     ω-vs-µ choice might be argued at first hand. **Highest-value acquisition.**
  2. **Alcouffe & O'Dell** (Hébert's ref. [36]) — the source he credits for the
     azimuthal `η_{p,q±1/2}` recursion (3.396)-(3.399); would say whether the
     azimuthal cell was ever intended to have edges in ω.
  3. **Reed & Lathrop** (Lathrop's ref. 7) — the cosine-picking lineage
     (Eqs. 27-29).

## The one-line answer

**For the cylinder's azimuthal closure the literature does NOT say** — no source
poses a weight in the azimuthal *angle*, and no source discusses the choice at
all: **BMC pose τ in the RADIAL DIRECTION COSINE µ (Eq. 74, per ξ-level,
declared "analogous to" the sphere's Eq. 42), Lathrop never treats the cylinder,
and Hébert prescribes the plain diamond τ ≡ ½ (Eq. 3.406) which is
variable-agnostic — but Hébert's own recommended cylindrical quadrature puts
every ordinate at the exact MIDPOINT IN ω of an equal-ω cell and weights it by
that cell's ω-measure, so on his rule τ ≡ ½ IS the barycentric weight taken in ω,
while BMC's cosine partition applied to that same rule drives τ out of the
admissible `[0,1]` and diverges with refinement.**
