---
name: lathrop-2000-angular-scheme-taxonomy
description: Lathrop 2000 NSE 134 is THE curvilinear angular-differencing taxonomy (5 schemes, sphere-only) and its Eq. 30 proves any tau != 1/2 is FIRST-order in angle; its Table II measures weighted diamond LOSING to plain diamond by up to 3.93x at grazing cosines
metadata:
  type: project
---

# Q68 — curvilinear SN angular differencing beyond the weighted diamond

Full deliverable: `scratch/q68_curvilinear_angular_differencing_survey.md`
(2026-08-12). This file is the negative-lookup row: what is settled, what is
still open, and where the sources are.

## The two spine sources — BOTH were already LOCAL and neither was in memory

**Lathrop, NSE 134(3):239–264 (2000), `10.13182/NSE00-A2114`.** LOCAL
(`Lathropp(2000)A Comparison of Angular Difference Schemes…`; printed p. = PDF
p. + 237). ⭐ **THE taxonomy.** Five schemes for the 1-D SPHERE, derived with
**no spatial differencing** (angular discretization → coupled ODEs in `r`,
Mathematica ODE solver, compared to a closed-form two-region absorbing analytic
solution). ⛔ `grep -ci cylind` = **0** — sphere only.

| scheme | eqs | unknowns/cell | angular order |
|---|---|---|---|
| step / constant discontinuous | (7),(8) | 1 | `O(Δμ)`; hybrid step+diamond tested → still 1st order on absorption |
| diamond (`δ=0`) | (18)+diamond; γ (21)/(22) | 1 | `O(Δμ²)` **iff node = midpoint** |
| weighted diamond (**Reed–Lathrop**) | (24)+(19)+(28) | 1 (ordinates are OUTPUTS) | `O(δΔμ + Δμ²)` |
| linear continuous (original S_N) | (9),(11),(12) | 1 | `O(Δμ²)` |
| linear discontinuous in ANGLE (Walters–Morel 1991) | (12)+(33), slope (32) | **2** | `O(Δμ²)` ptwise, **4th** on absorption/leakage |
| quadratic continuous (**new in this paper**) | (38)+(41) [or (43)]; (35)/(36)/(37) | **2** | perfect-derivative or `O(Δμ⁴)`; **4th** measured |

Cost: *"the higher-order schemes require at least twice as much arithmetic for
each Δμ"* (p.27).

**Larsen & Morel, "Advances in Discrete-Ordinates Methodology," Ch. 1 of
*Nuclear Computational Science: A Century in Review*, `10.1007/978-90-481-3411-3_1`.**
LOCAL (printed p. = PDF p. − 12). ⭐ **§1.5.1 is a 2-page review of exactly this
question.** ⚠ [[larsen-morel-2010-review-extraction]] already indexes it
(`pp.36-37 flux dip 1.92-1.93`) but carries **none of its verdict prose** — the
part that matters. "Indexed but not mined", not "missed".
`[SCAN p.49]`: *"Very few practical improvements beyond the elimination of the
flux dip have been made"*; angular DFEM *"has not happened"*; the LD-in-angle
scheme *"was found to be **less accurate** than the weighted-diamond scheme"*
because it abandons the accurate starting-direction flux; the **hybrid**
(quadratic-continuous first cell + LD after) wins. Their Eq. **(1.93)** is
Morel–Montry's τ verbatim. Ref **[26] = Dudziak, O'Dell & Alcouffe LA-7911-PR
(1979)** is confirmed as the **starting-direction-flux primary** (closes a
standing memory lead). Cylinder: only *"generalized to 2-D cylindrical geometry
[38]"*, i.e. M&M's own appendix.

## ⭐⭐ THE RULING on tau (answers the #235 floor and the 2.5–3.0× measurement)

**Lathrop Eq. (30) `[SCAN p.8]`** — third Taylor coefficient of the WD closure,
with `δ = 2τ − 1` (his Eq. 23/24 IS ORPHEUS's closure):

```
 −(1/4) Δμ [ 2δ(1 − μ_m²) + (1 − δ²) Δμ μ_m ]   ⟹  O(δΔμ + Δμ²)
```

*"only with `μ_m = μ̄` (`δ = 0`) is the truncation order `O(Δμ²)`."* Restated
twice in Sec. IV `[SCAN p.11]`, naming Reed's "optimum" WD explicitly.
⚠ Sidecar prints `−1/3`; **page says `−1/4`** (lessons L-010, re-confirmed).

⟹ **Every weighted diamond (Reed–Lathrop AND Morel–Montry) is FIRST-order in
angle. Plain diamond at the midpoint node is the only second-order option.** The
diffusion condition `c = Σ Δμ_m μ_m² = 2/3` (his Eqs. 29/54) is bought with
exactly that order: `δ=0` on an equal-interval mesh gives `c = 2/3 − Δμ²/6`
(Eq. 63), so the two properties are in **direct opposition on a fixed quadrature**.

**And it is MEASURED.** Table II `[SCAN p.14]`, max relative angular flux error,
test problem 1, ratio WD/diamond `[M]`: `+0.375 → 2.48×` (S8);
`+0.1875 → 2.14×`, **`+0.3125 → 3.93×`**, `+0.4375 → 1.86×` (S16). The losing
band is the **small-positive (outgoing near-grazing) cosines** and it **GROWS
with N** (2 of 8 → 5 of 16). ORPHEUS's `[M]` 2.5–3.0× sits inside this band.
⚠ Lathrop's own summary sentence says WD is better — true on the *average*
error; the direction-resolved max inverts it. **Read the table, not the abstract.**

**Clamp denominator now 6 of 6 primaries, NONE clamps τ.** Lathrop: *"`|δ| ≤ 1`
as long as `μ_m` is in `Δμ`"* — node-in-cell, i.e. `τ ∈ [0,1]`, nothing tighter.

## The (η,φ) 2-D closure — it EXISTS, and I nearly shipped a false negative

**Chaland & Samba, NSE 182(4):417–434 (2016), `10.13182/NSE15-38`. NOW LOCAL**
(I extracted it from the user's NSE archive + OCR'd it, 19 pp; printed p. = PDF
p. + 415). CEA/DAM. A **Cartesian product `(μ,φ)` angular mesh** — explicitly
*"instead of the standard nonconformal equi-solid-angle mesh"* — with **plain
diamond in BOTH** angular variables (Eq. 8) and **two** conservation recursions:
`β_{l+1} − β_l = −2μ̄Δμ` (μ) and `α_{m+1} − α_m = −cos φ̄ Δφ` (φ), both seeded at
0. 2-D angular sweep order from the angular characteristics (Eq. 7). ⭐ Their
parenthetical *"note that `β_l = 1 − μ_l²`"* is the cleanest published proof that
the conservation recursion hits the geometric value **exactly when the node is
the midpoint of an equispaced-μ cell**.
⛔ Disqualifiers: their own *"our scheme does not preserve the diffusion limit"*
(fixed only in App. B via Stone FE, order 2, diffusion limit deferred to a
follow-up paper); order 1 overall; mixed-coordinate frame motivated by ICF ray
effects, not 1-D accuracy; GMRES, no sweep. **Structural template, not a drop-in.**

## Actionable, beyond τ

**Lathrop Eq. (68)** `[OCR p.25]`: a midpoint-node equispaced rule integrates only
`P_0`/`P_1`, breaking anisotropic-source balance — the fix is to redefine
`P_k(μ_m)` as the **cell-averaged** `∫P_k dμ/Δμ = [(μP_k − P_{k−1})/((k+1)Δμ)]`,
which *"separate[s] the choice of quadratures for the streaming term and the
source term"*. Interpolating-polynomial variant credited to **Morel 1989 NSE
101:72 (LOCAL)**. ⟹ this is the missing third piece that makes "τ = ½ +
equispaced" a complete scheme. **Hébert (3.406) already ships τ ≡ ½** (plain
diamond in space AND angle for the cylinder) — the one openly-documentable
production code uses the scheme ORPHEUS measured as better.

## Negatives with their denominator (do not re-run)

- **No NAME for the (η,φ) azimuthal floor.** Named neighbours only: *flux dip*
  (M&M 1984), *ray effects* (Lathrop NSE 32:357, 1968), *thick-diffusion-limit
  failure* (spatial). Closest concept = **Hu & Azmy**, "regularity order w.r.t.
  the **azimuthal** angle caps angular convergence" — ANE 138:107199 (2020, **OA
  at `osti.gov/servlets/purl/1801122`**) + NSE 195:598–613 (2021,
  `10.1080/00295639.2020.1860634`) — but **2-D Cartesian, quadrature not stencil**.
- **No step / characteristic / exponential ANGULAR differencing** for curvilinear
  geometry anywhere searched.
- **Post-2010 is nearly empty on angular schemes.** Semantic-Scholar citation
  graphs of M&M 1984 (~30 citers) and Lathrop 2000 (**4** citers) fully classified;
  only Chaland–Samba qualifies. ⛔ Wang 2019 NSE 193:1339 *"Asymptotic Diffusion
  Limit of Numerical Schemes for the S_N Transport Equation"* is **SPATIAL**
  despite the BMC-like title. ⛔ Machorro 2007 JCP 223:67 DG on the 1-D sphere is
  **SPATIAL**. ⛔ Wu-Xie-Fischer 1999 NSE 133:350 (NOW LOCAL, OCR'd) is a
  **SPATIAL** nodal Green's-function method — angular part is the bare
  `β_{m+1} − β_m = −j W_m μ_m` recursion (j=1 cyl / j=2 sph).
- **NOT searched:** Scopus (403 all session), and **M&C / PHYSOR / ANS-Trans
  proceedings** — the real hole, since the LD-in-angle primary (Walters & Morel
  1991) is itself a proceedings paper.

## Bibliographic corrections `[DB]`

- **Bailey, Morel & Chang = NSE 165:149–169 (2010), DOI `10.13182/NSE08-66`** —
  *not* `NSE08-64` (that string floats around as a `research`-skill example).
  Worth checking `docs/refs.bib`.
- Walters & Morel 1991 has **two conflicting page pointers** in two local
  sources: Lathrop ref. 3 = "Vol. III, p. 13.2 3-1"; Larsen–Morel ref. [60] =
  "Sec. 11.1". Unresolved.

See [[morel-montry-tau-angular-cell-edges]],
[[reed-lathrop-1970-angular-truncation]],
[[larsen-morel-2010-review-extraction]],
[[user-nse-volume-archive]].
