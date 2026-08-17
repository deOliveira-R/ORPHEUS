---
name: morel-montry-tau-angular-cell-edges
description: Curvilinear SN angular closure — the 1-D cylindrical conservative streaming form (3 sources), the TWO cosines per azimuthal face, cylinder τ = BMC Eq.74 barycentric in the RADIAL cosine, diffusion-limit affine in η ALONE, no clamp, and the published typos.
metadata:
  type: project
---

Adjudication of ORPHEUS's Morel–Montry τ convention (2026-08-11, branch
`refactor/operator-strategy-layers`). Two deliverables, both with per-equation
scan verification:
- `scratch/q64_tau_edge_convention_literature.md` — the SPHERE / edge-convention pass.
- `scratch/q64_cylinder_closure_literature.md` — the **AZIMUTHAL (cylinder) ω-vs-η**
  pass; read it for anything cylinder-specific (§ "THE ω-vs-η RULING" below is its
  digest).

**Why:** an anisotropic cylindrical MMS put two τ conventions in conflict and could
not adjudicate (the τ ≡ ½ *control* beat every principled convention). The user
explicitly refused to let the MMS choose — "principled ≠ more accurate".

**How to apply:** this is the settled answer for any curvilinear angular-closure
question. Read the deliverable before re-deriving anything below.

## The ruling

**Cumulative-quadrature-weight cell edges, accumulated in the RADIAL direction
cosine, seeded at the level's inward endpoint — for BOTH geometries.** τ is the
ordinate's fractional position *in that radial cosine* inside its own cell.

| | edges | τ |
|---|---|---|
| sphere | **BMC Eq. (12)** `μ_{m+1/2}=μ_{m−1/2}+w_m`, `μ_{1/2}=−1`, `Σw=2` | **Eq. (42)** |
| cylinder (R-Z, per ξ-level) | **BMC Eq. (52)** same recursion, `μ_{1/2,n}=−√(1−ξ_n²)`, **weights renormalised per level to `Σ_m w̄ = 2√(1−ξ_n²)`** | **Eq. (74)**, declared "analogous to (42)" |

⟹ **The sphere/cylinder asymmetry is ORPHEUS's, NOT the literature's.** BMC derive
the cylinder as the sphere's analogue at every step (52↔12, 53↔15, 74↔42, 75↔41).
Independently confirmed for the sphere by **Lathrop 2000 p. 249** (`{μ_m, Δμ_m}`,
`μ_{m+1/2}=μ_{m−1/2}+Δμ_m`, `ΣΔμ_m = 2`).

**Refuted**: arithmetic mean of neighbouring ordinate cosines (nobody does this);
ω half-angle bisection (nobody does this either — see below).

## ⭐⭐ The adjudication the MMS cannot give — Lathrop 2000, two pages apart

- p. 245, after Eq. (30): *"only with `μ_m = μ̄` (δ = 0) is the truncation order
  `O(Δμ²)`"* — **τ = ½ has the BEST local truncation order.**
- p. 249→250, after Eq. (59): *"the standard angular diamond scheme, δ = 0, implies
  a cell-centered `μ_m = μ̄`, which … would not have `c = ⅔` and hence would not give
  the correct diffusion limit"* (`c ≡ Σ μ_m² Δμ_m`).

⟹ **The two criteria point opposite ways and the literature says so.** An MMS L2
error norm measures truncation; τ was chosen to zero BMC's β (Eq. 75), whose only
observable is the **flux dip / anomalous slope at the axis in a thick diffusive
problem**. The MMS is structurally blind to it. Lathrop's `δ` is τ affinely:
**τ = (1+δ)/2** (his Eqs. 23/24 = BMC 43/15).

**Cheap non-MMS oracle**: β is a pure function of the quadrature + the α/edge
recursions — **no solve needed**. Correct convention ⟹ β at round-off (BMC Table II:
`−4.12E−16` … `−1.55E−15`); wrong ones ⟹ `O(10⁻²–10⁰)`.

## The clamp — NO SOURCE HAS ONE

- Admissible range is `[0,1]`, and Lathrop says why: *"|δ| ≤ 1 as long as `μ_m` is
  in `Δμ`"* — a **consequence** of the ordinate lying in its cell, not an imposition.
- ⛔ **BMC's own worked S₂ gives `τ₁ = μ₁+1 = 1−1/√3 ≈ 0.4226 < ½`** (Eq. 47).
  ORPHEUS's `max(0.5, min(1.0, τ))` floor contradicts the source and re-introduces
  the β contamination τ exists to remove.
- The only negative-flux fixup in the corpus is **spatial** (Hébert 3.387–3.389).

## Hébert 2009 Ch.3 is a DIFFERENT scheme — τ ≡ ½, and he never defines edges

- Angular closure is the **plain diamond**, so labelled: **(3.406)** cylinder,
  **(3.431)** sphere. No τ, no weighted diamond, no Morel–Montry anywhere.
- ⚠ **The Eqs. ORPHEUS cites — (3.437)/(3.439) — are τ = ½ EXTRAPOLATIONS**
  (`φ_{n+1/2,i} = 2φ_{n,i} − φ_{n−1/2,i}` = BMC Eq. 16 at τ=½). They carry **zero**
  information about angular cell edges. They differ only in the *spatial* half
  (inward μ<0 vs outward μ>0 sweep).
- His cylinder defines only the **azimuthal** half-cosines `η_{p,q±1/2}` (3.398),
  and only because `α_{p,q±1/2} ≡ W_p η_{p,q±1/2}` (3.399). Attributed to
  **Alcouffe & O'Dell** (his ref. [36], NOT local).
- ⟹ the brief's "τ ≡ ½ control" is **not a control — it is Hébert's published scheme**.

## α vs τ — the premise of "do they share a partition?" is refuted twice

1. **α is NOT a geometric evaluation at a bisected ω.** Hébert derives (3.398) from
   *constant-flux preservation in an infinite non-absorbing medium*. The symbols
   `ω_{p,q±1/2}` appear only as labels for where the half-angle FLUX lives; no
   equation ever computes their values. Sole geometric anchor: (3.396), `η = 0` at
   the level's end point (direction lies in the ξ–μ plane).
2. **Lathrop states the relation explicitly** (Eqs. 22/25/26): `α_{m+1/2} =
   1 − μ_{m+1/2}²` **iff δ = 0**; in general `α = 1 − μ² + β` with
   `β_{m+1/2} − β_{m−1/2} = −δΔμ²/2`, `β_{1/2}=0`.

⟹ α (accumulates the **first** moment `Σwμ`) and the τ edges (accumulate the
**zeroth** moment `Σw`) are two different exact-integration conditions on ONE
partition. **Do not "fix" the apparent inconsistency by forcing them onto a common
geometric evaluation — that silently re-imposes τ = ½ and destroys α's conservation.**

## Two lineages — do not conflate

- **Reed & Lathrop**: choose the *cosines* to minimise truncation (Eqs. 27/28) ⟹
  fixes the quadrature; it then satisfies the diffusion condition `Σ Δμ μ² = ⅔`
  (Eq. 29) by construction, but integrates Legendre only through P₂.
- **Morel & Montry (1984) TTSP 13(5):615**: keep the standard quadrature (so
  anisotropic source moments stay right), **derive** τ from their own **Eq. (16b)**
  (sphere) / **(A2)** (cylinder) ⟹ **the quadrature is given and τ is derived —
  never free, never tuned.** ✅ **NOW LOCAL AND FULLY READ** (2026-08-11 third pass,
  §"MM-1984 PRIMARY" below) — it closes the clamp question: there is none.

## Notation traps

`μ`/`η` are **SWAPPED**: BMC & Hébert `μ` = ORPHEUS `η` (radial); their `η` =
ORPHEUS `ξ` (azimuthal); their `ξ` = ORPHEUS `mu_z`. Lathrop writes the **weight**
as `Δμ_m`. **The redistribution SIGN is not universal**: BMC (49) `+[…]/(rw)` ⟹
`α = α − wμ`; Hébert (3.400) `−(1/W)[…]` ⟹ `α = α + Wμ`. Check ORPHEUS's α sign
against ORPHEUS's own streaming sign, never against a quoted recursion.

## ⭐⭐ The RULE-AGNOSTIC predicate (what actually transfers) — 2026-08-11 follow-up

BMC Eq. (12)/(52) is **one solution**, not the condition. The condition, checkable on
any 1-D angular march `(η_m, w_m, e_{m±½}, τ_m)`:

| | predicate | source | order |
|---|---|---|---|
| **P0** | `Σw = \|range\|` and `Σ w η = 0` (makes α close) | Lathrop p.249; BMC 11/50 | — |
| **P1** | `c ≡ Σ w η² = ⅔` — *the diffusion condition* | Lathrop **Eq. 29/53/54** | **LEADING** |
| **P2** | `η_m = τ_m e_{m+½} + (1−τ_m) e_{m−½}` pointwise | **BMC Eq. 43 = Lathrop Eq. 23** | **FIRST** |
| **P3** | `η_m ∈ [e_{m−½}, e_{m+½}]` ⟺ `τ ∈ [0,1]` | Lathrop p.245 *"\|δ\|≤1 as long as μ_m is in Δμ"* | well-posedness |
| **P4** | `α_{m+½} = α_{m−½} − k w_m η_m`, seeded geometric | Hébert 3.397-3.399/3.424 | 1st moment |

**P2 in words = the ONE transferable sentence**: *τ is the barycentric coordinate of
the ordinate's marched cosine between its two cell-edge cosines* — equivalently, the
unique weight making the flux closure EXACT for ψ affine in that cosine. Silent about
edge provenance. **τ is determined, never chosen.** Its scalar shadow `β=0` is
strictly weaker `[M]`.

## ⭐⭐ THE β TENSION — two objects sharing a letter (resolution (i))

| | BMC Eq. 41/75 | Lathrop Eq. 25 |
|---|---|---|
| | `β = ½Σ η_m[α_{m+½}η_{m+½} − α_{m−½}η_{m−½}]` | `α_{m+½} = 1−η²_{m+½} + β_{m+½}` |
| type | ONE scalar (J⁽²⁾ contamination coeff.) | a SEQUENCE (α's pointwise defect) |
| zero iff | τ = Morel-Montry | **δ≡0 ⟺ τ≡½** (Eq. 26) |

**Mutually exclusive.** Mechanism, `[M]` **PROVEN by reproducing BMC Table I to every
printed digit**: BMC's `η_{m±½}` in β are the edges the SCHEME's τ IMPLIES via Eq. 43,
`ν_{1/2}=−1, ν_{m+½}=(η_m−(1−τ_m)ν_{m−½})/τ_m`. MM ⟹ ν ≡ cumulative-weight edges,
β=0. Diamond ⟹ ν drifts (S2: −0.1547 not 0; overshoots range 31 %), β≠0.
Lathrop's Eq. 57 oddness argument uses the TRUE edges, so β vanishes there — no
disagreement. **c=⅔ and β=0 live at DIFFERENT ORDERS** (BMC p.8: leading order is
accurate for ALL discretisations). Lathrop's "δ=0 ⟹ c≠⅔" is about the NODES δ=0
implies (midpoints); `[M]` midpoint c = 0.5/0.625/0.65625/0.664 at N=2/4/8/16.

⚠ Lathrop **Eq. 26's constant does not reproduce**: from his own Eq. 19+23 I get
`−δΔμ²`, page prints `−δΔμ²/2`. Structure (β ∝ δ) unaffected. L-010 class.

## ⛔ TWO REFUTATIONS of my own first-pass advice `[M]` 2026-08-11

1. **"(C) diverges from a missing weight renormalisation" — REFUTED.** With Eq. (52)
   applied correctly (`Σw̄ = 2 sinθ`), equispaced-ω nodes still give
   τ ∈ [−0.33,1.33] (M=8) → [−1.18,2.18] (16) → [−2.86,3.86] (32). **(C) violates P3**:
   `η_m = sinθ cos ω_m` is not inside the m-th *uniform-η* cell, and it worsens with
   refinement — matching the reported nan ladder. **Structural cause**: the arc cell's
   η-measure is `2 sinθ sin(ω_m) sin(Δω/2)` ∝ `sin ω_m`, NOT constant, while a
   trapezoid weight IS. ⟹ **BMC Eq. (52) is not a law — it is the statement that in
   THEIR quadrature the weight equals the cell's η-measure.** False by up to 5× for an
   equispaced-ω rule.
2. **"β is the oracle the MMS cannot be" — TRUE for the sphere, FALSE for a symmetric
   equispaced-ω level.** `[M]` β = round-off for (A),(B),(C) AND τ≡½ at M=8/16/32 —
   Mode-12 blind by Lathrop's own oddness. **Working gates instead**: (a) P3 as a
   PREDICATE (never clamp); (b) **ν-closure** — the τ-implied march must land on
   `+sinθ`; τ≡½ gives 1.0392/1.0097/1.0024 (M=8/16/32) while any derived τ closes
   exactly → the one cheap diagnostic separating diamond from a derived τ on this rule;
   (c) P2 pointwise residual; (d) P1 on the polar rule.

## ⭐⭐ THE ω-vs-η RULING — 2026-08-11 follow-up (deliverable `scratch/q64_cylinder_closure_literature.md`)

**No source poses the azimuthal closure in the ANGLE; no source discusses the
choice; and the strongest implicit evidence is nonetheless ω.**

- **Hébert's cylinder is §3.9.3 (printed 137-141), NOT §3.9.4** (= the SPHERE).
  Any brief citing "3.9.4 / 3.418-3.439" for the cylinder is mis-pointed.
- **(3.437)/(3.439) are NOT weighted** — they are (3.431) rearranged, and they are
  the SPHERE's polar pair. Cylinder counterparts = **(3.412)/(3.414)** (identical
  in the azimuthal half; only the *spatial* half flips). No weighted diamond
  exists anywhere in Hébert ch. 3.
- ⭐⭐ **`[M]` Hébert's OWN recommended cylindrical quadrature is equispaced-ω with
  UNIFORM weights.** (3.370) triangular and (3.375) product both reduce **exactly**
  to `ω_q = π(q−½)/M` (verified ≤4.4e-16, N=4…32); (3.371)/(3.376) `W_{p,q}=πW_p/M`
  is the cell's **ω-measure**. ⟹ his `τ≡½` (3.406) **IS** BMC-43/Lathrop-23 posed
  in ω — δ=0 *because the nodes really are the march-variable midpoints*. (An
  INFERENCE from quadrature+closure; he never says it.) Printed p. 133 recommends
  GL×Gauss-Chebyshev for 1-D cylinders; printed p. 135 = Eqs. (3.369)-(3.378).
- **BMC's cylinder τ is in the RADIAL COSINE µ, per ξ-level** — Eq. **(74)**
  (printed 160), declared *"analogous to … Eq. (42)"*. ω appears ONLY in the
  continuous Eq. (48) (`−(1/r)∂/∂ω[ηψ]`) and the coordinate line, then vanishes;
  the switch is never remarked on. β_cylindrical (75) is likewise all-µ.
- `[M]` **the brief's measure claim is CONFIRMED and it is what breaks the
  transplant**: an equal-ω cell's µ-measure is `2sinθ·sin(ω_q)·sin(Δω/2)` ∝ sin ω,
  spreading **5.0× / 10.2× / 20.4×** at M=8/16/32 (≈`2M/π`, WORSE with refinement)
  against a uniform weight. Two well-posed readings survive: **ω-barycentric ⟹
  τ ≡ ½ exactly**; **µ-barycentric on GEOMETRIC arc edges ⟹ `½+½cot(ω)tan(Δω/4)`
  ∈ [¼,¾]**. Sphere contrast: GL + Eq. 12 edges stays in `[0.39,0.61]`, terminus
  `+1.000000000000000`.
- **Lathrop 2000 has ZERO cylinder content** ("cylind" count = 0) — sphere-only;
  cylinder use is by analogy through BMC §V.
- **Lathrop treats τ≡½ as a DEFICIENCY, twice** (settles the brief's "which
  direction does his argument run"): printed **245** — the diamond gives
  `ψ_m = ψ(r,µ̄)`, *not* `ψ(r,µ_m)`, error **O(δΔµ)**, remedy = Eq. (24); printed
  **249-250** — δ=0 *implies* midpoint nodes, and **Eq. (63)**
  `c = ⅔(1−1/N²) = ⅔ − Δµ²/6` proves they miss the diffusion condition. ⟹ Eq. 63
  is the CLOSED FORM of the previously-measured midpoint-c ladder
  (0.5/0.625/0.65625/0.6640625 at N=2/4/8/16) — no re-measurement needed.
- **CLAMP: still nobody.** Range is `[0,1]` in all three (BMC p.152 & p.157 *"any
  value between zero and one"*; Lathrop `|δ|≤1` as a CONSEQUENCE of node∈cell).
  Positivity is analysed only for the **SPATIAL** diamond (Hébert 3.387-3.389,
  `Δx Σ < 2|µ|`); Lathrop's sole negative-flux mention is explicitly *spatial*
  (printed 262). **The `τ≥½ ⟺ |(1−τ)/τ|≤1` stability argument is in NO source** —
  BMC even print the amplifying recurrence (54) and say nothing. `[M]` BMC Eq. (47)
  reproduced: S₂ GL `τ₁ = 1−1/√3 = 0.42265 < ½`, `τ₂ = 0.57735`; exactly HALF the
  set is sub-½ at every order.
- ⛔ **Three PUBLISHED typos, all verified on the rendered page** (L-010): BMC
  Eq. (50) `α_{m+1/2,n}=α_{m+1/2,n}−µw` (RHS ⟹ `m−1/2`), Eq. (52) same slip for µ,
  and printed p. 156's coordinate line `η = sinθ cos ω` (⟹ `sin ω`, else η≡µ).

## Cylinder specifics `[M]` + `[VERIFIED-ON-SCAN p.72]`

- **α's target = `w_gl × (tangential cosine at the boundary)`** — Hébert (3.399)
  verbatim. Exactness = (3.398): `w_m η_m = w_gl(ξ_{m+½} − ξ_{m−½})`.
- **ORPHEUS's T3 α IS Hébert's α exactly.** Running (3.398) on equispaced-ω +
  trapezoid weights gives `ξ̃_{m+½} = κ·ξ(ω_{m+½})`, `κ = Δω/(2sin(Δω/2))`. So κ is
  NOT an error in α — it is the ratio between the **recursion-defined edge (which is
  what α *is*, in both sources)** and the geometric arc edge (**which no source uses**).
- **`κ−1` IS the cylinder's Lathrop-β** (α's departure from its exact target) —
  CONFIRMED structurally — **which is exactly why driving κ→1 is wrong**: it drives
  δ→0, i.e. τ→½.
- **Neither source faces the two-weights conflict**: they let the **EDGE float** and
  keep the weight. Demanding α-exactness at geometric edges forces
  `w_m = 2w_gl sin(Δω/2)` (nodes untouched, uniform `1/κ`) — confirmed — but that
  breaks `Σw = π` / P0. **Do not rescale the weights.**
- `[M]` **(A) IS (B) in disguise**: (A)'s interior edges = `cos(Δω/2) × (B)`'s arc
  edges to 1e−16, end cells stretched. ⟹ prefer **(B)**, which satisfies P3 at all M
  and has the closed form **`τ_m = ½ + ½ cot(ω_m) tan(Δω/4)`** (verified 6.7e−16).
- `[M]` **At S₈ Gauss–Legendre FOUR of eight MM τ are below ½**
  (0.3923, 0.4591, 0.4809, 0.4942) — the clamp corrupts half the set, not one value.

## ⭐⭐ 2026-08-11 SECOND PASS (`scratch/q64_attempt2_lit_check.md`) — the CYLINDER is now fully sourced

### The 1-D cylindrical CONSERVATIVE STREAMING form — SETTLED, three independent sources

`Ω·∇ψ = (η/r)∂(rψ)/∂r − (1/r)∂(ξψ)/∂ω` (ORPHEUS letters) is published verbatim:
**Hébert (3.157)** p. 92 `[SCAN]` · **Bell & Glasstone 1970** p. 58 display + **Table 1.2**
"infinite cylinder, axial symmetry" row p. 59 `[SCAN]` · **BMC Eq. (48)** p. 156 `[SCAN]`.
Non-conservative twin = **Hébert (3.155)** / B&G Table 1.2. Consistency: Hébert states
`∂η/∂ω = μ` (ORPHEUS `∂ξ/∂ω = η`), which makes the two forms identical, not approximate.
`∫_{4π} ∂(ηφ)/∂ω = 0` = **Hébert (3.158)** = B&G's "removes the last two terms".
**ω DECREASES along the path** (Hébert 3.149, with his one-clause reason), sweep starts
at ω=π (3.407 + 3.373 + Fig. 3.27), Carlson–Lathrop order "μ_m increases on [−1,1]"
⟹ **ORPHEUS's η-ascending / ω-descending storage IS the literature's traversal.**
⛔ `SILENT`: Pomraning 1989 (declines the cylinder BY NAME — never cite it for this),
Ligou Ch.8, LA-3186/4058, Lathrop 2000, Sanchez–Ganapol 1983.

### ⭐⭐ TWO cosines live at every azimuthal face — this reconciles the apparent contradiction

| at face `m±½`, level `n` | defined by | used for |
|---|---|---|
| **azimuthal** cosine (ORPHEUS `ξ`) = `α/W` | conservation recursion **BMC (50)** / **Hébert (3.398)** / **LA-3251 (3-10)** | the **streaming coefficient** |
| **radial** cosine (ORPHEUS `η`) | cumulative-weight recursion **BMC (52)** | **τ** (Eq. 74) and **β** (Eq. 75) |

⟹ "the face coefficient is ξ" and "τ is barycentric in η" are BOTH true, about
DIFFERENT numbers at the SAME face. **Stop treating it as a conflict.**

### ⭐ THE CYLINDER τ IS **BMC Eq. (74)**, not 42/43 — and it is barycentric in the RADIAL cosine

`τ_{m,n} = (μ_{m,n} − μ_{m−½,n})/(μ_{m+½,n} − μ_{m−½,n})`, printed p. 160 `[SCAN]`,
BMC's `μ` = **radial** cosine. ⟹ **never in ω.** ORPHEUS's code citing "BMC 43 =
Lathrop 23" is citing the **SPHERE** pair — same content, wrong geometry's number.
**BMC Eq. (74) is the LONE cylinder τ in the corpus** (Hébert has no τ at all).

### ⭐⭐ THE CRUX (2(c)) — the cylinder closure is graded on ψ affine in the RADIAL cosine ALONE

Three printed places, all `[SCAN]` pp. 150/158/159:
1. **BMC (61)+(62)**: `Ω_{m,n} = μ e_r + ξ e_z`, `∇ = ∂_r e_r + ∂_z e_z` — **NO azimuthal
   component**. ⟹ `ψ^{(1)}` affine in (radial, axial); in 1-D, **radial alone**.
2. **BMC Eq. (1)**: `ψ_m = φ/4π + (3/4π) J_r μ_m`, `J_r` glossed *"radial component of
   the current"* — one first-order mode, radial.
3. **BMC (69)** is a PAIR of conditions (radial-edge, axial-edge); **(70)** kills the
   axial one identically *"because `ξ_{m−½,n} = ξ_{m+½,n} = ξ_{m,n}`"* (no cell
   structure on a level). **There is NO azimuthal condition anywhere in the paper.**

⟹ **`REFUTED`: no source treats the azimuthal cosine as an independent first-order
mode**, and two state the symmetry forbidding it — **Hébert (3.153)** `φ(ρ,ξ,ω) =
φ(ρ,ξ,−ω)` `[SCAN]` and **LA-3251 printed p. 31** `[SCAN]`: *"In one-dimensional
cylindrical geometry the flux is symmetric in ξ and η, so that only a quadrant,
containing the entire range of μ, is needed"* (his cyl letters: η = azimuthal). Hence
`J_φ ≡ 0`. **Grade a cylindrical closure on "affine in η", never on "affine in (η,ξ)".**

### The τ EXACTNESS claim is TWO statements that are the same property

(i) pointwise: under **(43)** — *"Eq. (15) will exactly relate the cell-edge and
cell-center fluxes when the angular flux assumes the linear form defined by Eq. (1)"*;
(ii) global: after **(74)** — *"first-order consistency in the diffusion limit"* =
*"preserve the Galerkin diffusion approximation"* (75), and MM is *"the only method …
that forces the β factor to be zero for any standard quadrature set"* (p. 161).
They coincide because **(35)/(61)** make the diffusion-limit `ψ^{(1)}` affine in the
radial cosine. ⚠ **τ is irrelevant at LEADING order** — BMC p. 159 after (66):
*"true for any choice of weighting factors."*

### The CLAMP — `REFUTED` three ways `[SCAN]`

1. BMC under **(53)** and **(15)**: *"the weighting factors can take on **any value
   between zero and one**"* — `½` is an interior landmark (diamond), not a floor.
2. **Eq. (47)** prints `τ_1 = μ_1 + 1`, `τ_2 = μ_2`; **Fig. 1** prints `μ_1 = −(1/3)^{1/2}`
   ⟹ **`τ_1 = 1 − 1/√3 ≈ 0.42265 < ½`**. Both halves are on p. 155, so the number is
   fully sourced (formula + ordinates), though not printed as a number.
3. After (47): `τ ≡ ½` *"would force the quadrature set to be the midpoint rule …
   However, the midpoint rule is not a good choice of quadrature sets in general."*

### ⭐ Hébert's OWN cylinder quadrature IS ORPHEUS's rule — (3.369)–(3.373) `[SCAN p.135]`

`ω_{p,q}` = **cell-centred equispaced**, `Δω = π/(N−2p+2)`, ω INCREASES with q;
`W_{p,q} = πW_p/(N−2p+2)` ⟹ **all within-level weights EQUAL**; `ΣΣW = π`; zero-weight
starting points at **ω = π** (3.373, TOP half-index). ⟹ `W_{p,q}/W_p = Δω`, so (3.398)
is `η_{q+½} = η_{q−½} + Δω·sinθ·cos ω_q`.
`[M]` its output = **`κ ×` geometric arc value at EVERY interior face**,
`κ = Δω/(2 sin(Δω/2))` (M=8: `1.0064545428`, agreement 1e-12), closing at `1.1e-16`.
⟹ **κ is an `O(Δω²)` conservative-vs-geometric gap, not an error; κ→1 breaks
conservation.** And Hébert pairs this quadrature with the **plain diamond (3.406)**.
⟹ Hébert-scheme and BMC-scheme are two complete DIFFERENT schemes; ORPHEUS mixes
Hébert's quadrature with BMC's τ — legitimate only through **(43)/(74)**'s
edge-agnostic form, **NEVER through (52)** (whose hidden premise, weight = radial-cosine
cell measure, is false by up to ≈5× here — the real cause of the τ∉[0,1] ladder).

### ⚠⚠ THE τ SIGN TRAP — the closed form flips with the MARCH ORIENTATION `[M]` 1e-16

```
index marches WITH ω      (Hébert's q)                 → τ = ½ − ½ cot(ω_m) tan(Δω/4)
index marches WITH the radial cosine = AGAINST ω
        (BMC Eq. 52 seed −√(1−ξ²) → +√(1−ξ²) = ORPHEUS) → τ = ½ + ½ cot(ω_m) tan(Δω/4)
```
Same SET, reversed ORDER (M=8: 0.2524…0.7476), both in [0,1] at all M.
⟹ **ORPHEUS's `+` is correct for its η-ascending storage.** A symmetric fixture cannot
see the flip.

### New PUBLISHED typos found (L-010 family) — all `[SCAN]`, OCR faithful

- ⛔ **BMC p. 156, under Eq. (48): *"η = sin θ cos ω"*** — must be `sin ω`. Refuted 4
  ways (μ≠η; `μ²+η²+ξ²≠1`; `∂η/∂ω=μ` fails; Fig. 2). **The most dangerous trap in the
  corpus**: taken literally it makes the azimuthal-face coefficient the radial cosine.
- ⛔ **BMC Eq. (50)** printed self-referentially (`α_{m+½}=α_{m+½}−μw`); RHS = `m−½`.
  Same family as (52). Sphere twins (11)/(12) are printed correctly.
- ⛔ **Stacey (9.150)** p. 344: `sinθ/4` must be `sinφ/r` — **dimensionally impossible
  as printed**. The only local source that *reads* like it contradicts the claim.

### Provenance closed / still open

- ✅ **Morel & Montry 1984 VERIFIED via CrossRef**: TTSP **13**(5), **615–633**, DOI
  **`10.1080/00411458408211661`**, cited-by 19. Prior notes were right. ~~Still not
  local~~ ⟹ **NOW LOCAL, OCR'd, and READ END-TO-END** (2026-08-11) — see the
  MM-1984 PRIMARY section below; deliverable `scratch/q64_morel_montry_findings.md`.
- ⛔ **Alcouffe & O'Dell (Hébert ref. [36], named at his p. 138 as the author of the
  cylinder `η_{p,q±½}` construction) — UNRESOLVED after 7 queries across 4 databases.**
  Hébert's Ch.-3 reference list is NOT in the local extract. ⟹ **the primary source for
  ORPHEUS's cylinder angular-cell-edge construction has never been read.** Do not cite
  a guessed venue. See [[lessons-L012]].
- **Zotero DOWN** (port 23119 refused) — no annotations checked; re-query before any
  theory page lands.

## ⭐⭐⭐ MM-1984 PRIMARY — read end-to-end 2026-08-11 (`scratch/q64_morel_montry_findings.md`)

**Morel & Montry 1984 TTSP 13(5):615-633 is LOCAL + OCR'd. Printed p. = PDF p. + 613.**
19 printed pp., **zero tables**, eqs (1)–(19) + (1a-d),(4a-b),(6a),(7a),(16a-b),(17a-b)
+ (A1)–(A4). `[SCAN]` verified: PDF 4 / 11 / 12 / 17 / 20.

### ⭐⭐ THE ENDPOINT RULING — the "BMC-43 vs R&L-16 conflict" is a DECLARED TRADE, not an open question
**PDF p. 4 / printed 617 `[SCAN]`**: *unlike the Reed–Lathrop scheme ours works with
**any** S_N quadrature set; however, like the standard diamond scheme, **"our scheme is
only first-order accurate."*** M&M **cite R&L as their ref. 1** and reject them at
printed 616 (R&L's sets can't integrate degree > 3 ⟹ breaks conservation for
anisotropic scattering). ⟹ **barycentric-in-μ WINS; second-order-in-angle is the
price the primary knowingly pays.** Nothing licenses clamping τ toward ½.

### The clamp denominator is now **5 of 5 primaries, NONE clamps** (M&M · R&L · BMC · Hébert · Lathrop)
τ ∈ (0,1) automatically (monotone edge march + node-in-cell); no *tighter* interval
is ever proposed, in either geometry. `[½,1]` stays Grant's, on the SPATIAL weight.

### ⭐⭐ THE CYLINDER — MM's appendix is a REAL derivation, and it RATIFIES the 6.4 partition
Printed **632** in as many words: in the cylinder the cell-edge **cosines** are *not*
weight-partial-sums as in the sphere — the cell-edge **AZIMUTHS** are; the edge cosines
follow from cell-centre polar + cell-edge azimuth. Eqs. `[SCAN p.20]`:
`(A1)` `ψ=(1−τ)ψ_{m−½}+τψ_{m+½}` · `(A2)` τ barycentric in the **radial cosine** ·
`(A3)` `μ_{m+½}=sinθ_ℓ cos φ_{m+½}` · `(A4)` `φ_{m+½}=φ_{m−½}+πC_ℓW_{ℓ,m}`, `φ_{1/2}=−π`.
⟹ **partition in ω by weight; barycentre in μ.** Equal weights ⟹ (A4) IS equispaced-ω
midpoints ⟹ `[M]` `τ_m = ½ − ½cot(φ_m)tan(Δφ/4)` (machine-exact vs A2–A4, M=2…64),
**independent of sinθ_ℓ**. ORPHEUS's `+` = the same set, march reversed.
`[M]` `τ₁→¼`, `τ_M→¾`, `(τ₁−½)/w = −M/4` exactly; end/mid μ-width ratio ∝ 1/M².
⟹ **ORPHEUS's observed ¼ IS M&M's own formula. They neither flag nor bound it.**
⚠ **The appendix has NO cylinder β / diffusion-limit analysis** — recipe + "we tested
it, the dip is similarly eliminated". Every β claim in the paper is SPHERICAL.

### ⛔ BMC line 657 ("forcing β to zero DETERMINES the MM weights") is FALSE as stated
M&M's defining condition is **per-ordinate edge COINCIDENCE**: make the scheme-implied
edges (4a) equal the standard weight-defined edges **(15)** `μ_{m+½}=μ_{m−½}+2W_m`,
`μ_{1/2}=−1` (their ref. 3 = Bell & Glasstone) — **N equations**, inverted ordinate-by-
ordinate ⟹ **(16b)**. β=0 is the **COROLLARY**, proved by parity in (17a)–(19) (δ_m odd,
`α_{m+½}+α_{m−½}` even). `[M]` β is ONE scalar: root-finding τ₁ on Gauss-S8 with the
other seven randomised gives β=O(1e-16) at `‖τ−τ_MM‖_∞ = 0.238 / 0.308 / 0.242`.
⟹ **β=0 is NECESSARY, not sufficient — an (N−1)-dim family.** See [[lessons-L013]].

### ⭐⭐ THE SEED — unit gain, no damping, sign-determines the dip
`(1a)` `[SCAN]` `α_{m+½}=α_{m−½}−μ_mW_m`, **both** `α_{1/2}=0` **and** `α_{N+1/2}=0`
printed; **NOTHING on march direction, conditioning, or error growth.** But `(14)`
`[SCAN p.11]` admits an arbitrary edge-march seed `μ_{1/2}=μ_s`, and for Gauss-S₂ +
diamond β goes negative for `μ_s < −2/√3`. `[M]` reproduced exactly: **`β(μ_s) = μ_s +
2/√3`** — linear, unit gain. Their instrument is the *effective starting cosine*
`(10)` `ψ̃_s = φ̃ + 3J̃μ_s` (Figs. 4/6/8; healthy = μ_s → −1 in the interior; the
S₂-diamond failure bottoms at ≈ −1.35). Underestimate ⟹ β<0 ⟹ **dip**; overestimate ⟹
β>0 ⟹ no dip but the diffusion limit is lost. Printed 629: spatial truncation is
`Δr²/r²` for weighted fluxes vs `Δr²` for the starting flux ⟹ coarse-mesh imbalance.

### `[M]` The endpoint defect is ALREADY in M&M's own SPHERE (Gauss + Eq. 15)
τ₁ = 0.4226 (S2) → 0.3923 (S8) → 0.38971 (S64): **saturates, never → ½**; `(τ₁−½)/(2W₁)`
= −0.077 / −1.06 / −61.85 ⟹ **O(N) divergence, i.e. R&L's degradation, in the primary's
default configuration.** Interior → ½ (0.500099 mid-level at S64). Their other claims
reproduce: β_MM = 0 at N=2…32; β_step = 0.577→0.043 (S2→S32, "converges to zero");
β_diamond = +0.1547 at S2 (their Fig. 3 "β is positive"), −2.7e-3 at S4.

### POSITIVITY — no authority for a half-angle gate
No statement anywhere on `ψ_{m±½} ≥ 0`. The only negativity analysed is of **D** (7a)
via β<0: `r_D = −2β/σ_t` is where D blows up and the slope reverses — **that IS the
dip**, a scalar-flux local minimum, not a negative flux. (Consistent with the earlier
finding that all corpus positivity work is SPATIAL.)

### Errata (L-010) — 2 print, 2 OCR
⛔ **Eq. (1) is printed with the two ψ half-indices CROSSED** (`α_{m+½}ψ_{m−½} −
α_{m−½}ψ_{m+½}`); their own (4) and (6a) are uncrossed, and only uncrossed telescopes
under the zeroth moment (⟹ Eq. 5 = exact continuity). **PRINT typo, scan-verified.**
· `[OCR]` (12)/(13): sidecar reads `dφ/dr`, page reads `dψ/dr`; `[M]` algebra decisive
— with `dψ_s/dr`, (11)+(8) satisfies (12) and (11)+(9) gives `Q(1+2β)`; with `dφ/dr`,
(11)+(8) gives 0 ≠ Q. · ⚠ printed 630 *"ξ = sin θ"* (ξ is conventionally the AXIAL
cosine; inert for A2–A4, but never import M&M's ξ into a B&G/BMC context). · printed
630 Summary *"relative to the starting fluxes"* ⟹ "weighted". · `[OCR]` cover page
"J. E. Horel" ⟹ Morel.

### Notation map (M&M → ORPHEUS)
sphere `μ` = streaming cosine · cylinder `μ = sinθ cos φ` = **radial** = ORPHEUS `mu_x`
· `φ` = azimuth = ORPHEUS `ω` · `θ_ℓ` = level polar · `ξ`-level = θ-level.
`Σ W_m = 1` (NOT 2) ⟹ **weight = HALF the cell μ-measure**, hence (15)'s `2W_m` and
(1)'s `2/r` prefactor. ⚠ M&M's letters ≠ BMC's (BMC radial = η, azimuth = ω).

### ⭐ Lead opened: the "Miller–Alcouffe" attribution mismatch
In-text credits *"Miller and Alcouffe"*; **ref. 2 is D. J. Dudziak, R. D. O'Dell &
R. E. Alcouffe, "Transport and Reactor Theory," LANL LA-7911-PR (July 1979)** — no
Miller. ⟹ candidate for the unresolved **Alcouffe & O'Dell** primary (Hébert ref. 36).
`LA-####-PR` T-1 progress reports ARE full-text on OSTI. ⛔ **UNREAD/UNVERIFIED** —
OSTI client hit `ConnectionResetError` twice this session; a web search conflated
quarters. Resolve the exact record before citing. See [[lessons-L012]].

See [[sphere-sn-pole-closure-canonical]], [[curvilinear-sweep-directness-ruling]],
[[reed-lathrop-1970-angular-truncation]], [[lessons-L010]], [[lessons-L012]],
[[lessons-L013]], [[lessons-L014]].
