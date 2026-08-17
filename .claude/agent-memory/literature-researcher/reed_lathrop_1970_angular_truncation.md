---
name: reed-lathrop-1970-angular-truncation
description: Reed & Lathrop 1970 NSE 41(2):237-248 — the ANGULAR truncation-error analysis; what "R&L weighted diamond" actually is (ordinates are OUTPUTS), the tau=1/2+O(w) second-order precondition, the two solve-free gates, and the cylinder's separated weight-vs-Delta-mu form
metadata:
  type: reference
---

# Reed & Lathrop (1970) — angular difference scheme, truncation error

W. H. Reed and K. D. Lathrop, *NSE* **41**(2), 237–248 (1970),
DOI `10.13182/NSE70-A20710` — ✅ CrossRef-verified 2026-08-11.
LOCAL: `scratch/literature/07-Truncation Error Analysis…pdf` (13 pp.);
sidecar same stem. **printed = PDF + 235.** Full extraction:
`scratch/q64_reed_lathrop_findings.md`.
⚠ Sidecar OCR misreads half-integer subscript `1` as `3` on PDF p. 3
(Eqs. 2/4a/4b) — the PRINT is right; verify on the render.
Table I body is in `.mocr.json` `pages[11].tables[0]` (L-008).

## ⭐ The ruling: "R&L weighted diamond" ≠ "Morel–Montry weighted diamond"

Both use the same closure `N = τ N_{m+1/2} + (1−τ) N_{m−1/2}` (**Eq. 4b**,
printed p. 238) and the same barycentric node `μ_m = τ μ_{m+1/2} + (1−τ) μ_{m−1/2}`
(**Eq. 3b** = BMC 42/43 = predicate P2), both credited to **Grant 1968**.
The fork is one extra equation — the angular **consistency condition**:

| | R&L 1970 | Morel–Montry / BMC (ORPHEUS) |
|---|---|---|
| (13a) `α_{m+1/2} − α_{m−1/2} = −μ_m w_m` | yes | yes |
| **(13b) `(1−τ)α_{m+1/2} + τα_{m−1/2} = (1−μ_m²)/2`** | **yes** | **no** |
| (13c) barycentric τ | yes | yes |
| ordinates `μ_m` | **OUTPUTS** (Eq. 14 quadratic, marched from `α_{1/2}=0`) | inputs |

**"Fixes a quadrature" (Lathrop 2000's phrase) decoded:** weights in, ordinates
out. Eqs. **(14a)/(14b)** printed p. 240 are quadratics in `μ_m`. Consequences the
authors state themselves: `Σ w μ² = ⅓ Σ w` holds automatically for ANY symmetric
interval set (property 2, needs `α_{1/2}=α_{M+1/2}=0`); the abstract says exact
"to second order in the direction cosines"; **Conclusions p. 248: 4th and higher
even moments NOT correctly integrated, "<3 % for the fourth moment and S₈ Gauss
weights."** That last sentence IS Lathrop 2000's objection, in the primary.

## The two solve-free gates this paper hands you (the actionable output)

1. **`(τ_m − ½)/w_m` must stay BOUNDED under N-refinement.** From property 3 +
   Eqs. (15)/(16) printed p. 241: the `w²` coefficient
   `[(1−τ)²α_{m+1/2} − τ²α_{m−1/2}]/w_m` is `O(1/w)` — angular scheme degrades
   2nd→1st order — unless the node is the cell midpoint to `O(w²)`, i.e.
   `τ = ½ + O(w)`. ⚠ ASYMPTOTIC: forbids `τ−½` not shrinking like `w`, not a
   fixed-N value.
2. **The (13b) residual must be `O(w²)`.** `R_m ≡ (1−τ)α_{m+1/2} + τα_{m−1/2}
   − (1−μ_m²)/2`; a non-zero `R_m` leaves the angular streaming coefficient wrong
   by `2R_m/r` multiplying `∂ψ/∂μ` (from Eq. 7 with Eq. 12). Calibration: R&L say
   plain diamond's residual is `O(Δμ²)` (p. 246). **Pointwise ⟹ NOT blind on a
   symmetric level**, unlike β / ν-closure.

## `[M]` My reproduction (2026-08-11) — reproduces properties 1 & 2 to 1e-14

R&L's own determined τ ∈ `[0.4397,0.5603]` at S8-Gauss → `[0.4827,0.5173]` at
S128. Per-step `|(1−τ)/τ| ≤ 1.27`; transient peak (running product) **1.5–1.8**,
always at the half-level; **end-to-end product EXACTLY 1**.
⭐ That last is a THEOREM not a fluke: property 1 forces `τ_{M+1−m} = 1 − τ_m`, so
factors pair to 1. **Holds for ANY barycentric τ on a symmetric partition with
symmetric ordinates — including ORPHEUS's.** Cheap invariant to assert.
⚠ A ONE-SIDED α march across the whole level is ill-conditioned (S64 Gauss: loses
symmetry, gain → 1e9). R&L prescribe TWO marches, each seeded at its own `α=0`
end, meeting at μ=0 — they frame it as symmetry (p. 240), but conditioning follows.

## Answers to the six 2026-08-11 questions (q64 amplification/seed brief)

- **Truncation error separated spatial vs angular?** YES (Eq. 7's two braces).
  Spatial with Eq. (11) `a = ½ + (1/6)Δr/(r_{i+1}+r_i)`: `(Δr²/r)∂²ψ/∂r²` vs plain
  diamond's `(Δr²/r²)∂ψ/∂r`.
- **`(1−τ)/τ` amplification / recurrence stability?** ⛔ **NOT ADDRESSED.** No
  `amplif|accumul|stabil|oscillat|monoton` anywhere. τ enters only as a LOCAL
  truncation coefficient, never as a recurrence gain.
- **Positivity / negative flux?** ⛔ **NOT ADDRESSED.** Sole mention: they REMOVED
  DTF-IV's negative-flux fixups so the comparison measured differencing (p. 246).
  `τ ∈ (0,1)` (p. 238) is node-containment (P3), open interval, no clamp.
- **Starting direction?** PARTIAL. `α_{1/2}=0` seeds; **Eq. (17) `μ_1 = −1 + w_1/2
  + O(w_1²)`**; Eqs. (18)/(19) are an INDUCTION carrying the near-midpoint property
  from the seed through every ordinate, hinging on a root-sign rule (negative root
  for μ<0, positive for μ>0). **But it propagates the ORDINATE property, not a
  FLUX error** — the μ=−1 seed flux equation is never written or analysed.
- **Geometry:** §II.A 1-D spheres, §II.B **3-D CYLINDERS** (pp. 242-243, Eqs.
  20-28), §II.C 3-D spheres. Cartesian explicitly excluded. **Numerics 1-D sphere
  ONLY.** ⟹ contra Lathrop 2000, this paper DOES carry the cylinder.

## ⭐ Cylinder: weight `w` and cell width `Δμ` are kept SEPARATE

Eq. **(28b)** printed p. 243: `[(1−τ)α_{m+1/2} + τα_{m−1/2}] Δμ = η² w`, with
`η² = 1 − μ_m² − ξ_l²`; α recursion (28a) uses `w`. The sphere collapses them only
because p. 238 defines `w_m` as "**normally** `μ_{m+1/2} − μ_{m−1/2}`" (note the
hedge). **Primary-source licence for weight ≠ cell-measure in a cylinder** — the
crux behind "BMC Eq. 52 is NOT a law". Also: seeds `α_{1/2}=α_{M+1/2}=0` **per ξ
level**, μ-interval count may differ per level, ξ levels otherwise FREE.
Cylinder angular derivative is in the RADIAL cosine: Eq. (21)
`… + (η²/r)∂ψ/∂μ`, via `∂ψ′/∂ω = −η ∂ψ/∂μ`.

## ⛔ NOTATION: R&L's cylinder cosines are a 3-CYCLE off BMC/Hébert

| axis | ORPHEUS | **Reed & Lathrop** | BMC / Hébert |
|---|---|---|---|
| radial | `mu_x` | **μ** | η |
| azimuthal | `mu_y` | **η** | ξ |
| axial | `mu_z` | **ξ** | μ |

Every symbol is occupied by a different axis. Transcribing Eq. (28)
symbol-for-symbol against a BMC-conditioned reading gives the wrong equation.
Also: R&L's `α′ = (A_{i+1} − A_i)α` is AREA-SCALED — not ORPHEUS's `alpha`.
Their `τ`, `α`, `w` otherwise map 1:1 (their α recursion = Hébert 3.398, same sign).

## Table I (printed p. 247) — the δ=0 ⟹ midpoint-nodes experiment, run in 1970

6-group ANL-7416 spherical benchmark, `k_eff = 0.99597`, 40 intervals. % error:

| S_N | Gauss, diamond | **Gauss WEIGHTS only, node = interval average** | Gauss wts + Eq.14 μ, diamond | + angle wt | equal wts + Eq.14 μ, both wts |
|---|---|---|---|---|---|
| 2 | 3.577 | **13.87** | 3.577 | 2.094 | 2.096 |
| 4 | 0.921 | **4.29** | 0.965 | 0.737 | 0.407 |
| 8 | 0.259 | **1.20** | 0.225 | 0.201 | 0.097 |
| 16 | 0.070 | — | 0.055 | 0.053 | 0.024 |
| 32 | 0.019 | — | 0.014 | 0.014 | 0.006 |

Col 2 = keep Gauss weights, move each node to its interval centre so τ≡½ is exact;
paper's verdict: "does not satisfy the diffusion condition and leads to the poor
results shown" (p. 246). **4.6× worse at S₈.** Effect decomposition at S₈:
quadrature 2nd moment **5.3×**; closure weight τ only **11 %**; spatial weight
**nil** on `k_eff`. ⟹ in R&L's configuration accuracy is dominated by the
QUADRATURE, not by τ and not by any recurrence gain — **but their τ never leaves
`½±0.06`, so their configuration CANNOT exhibit a 9× transient; silence ≠ absence.**

⚠ **Flux-dip attribution:** in THIS paper the central flux dip is cured by the
**SPATIAL** weight (Eq. 11), which has "very little effect on the eigenvalue";
the ANGULAR weight is what improves `k_eff` (p. 246). Whether Morel–Montry's dip
is the same object is NOT decidable here.

## Clamp denominator (extends L-012 discipline)

**4 of 4 sources that define τ — BMC 2010, Lathrop 2000, Hébert 2009, Reed &
Lathrop 1970 — clamp NOTHING.** The `[½,1]` interval in this lineage is **Grant's,
on the SPATIAL weight `a`** (R&L footnote 8, printed p. 239: Grant's sign-dependent
weight is needed "only to keep `a` between ½ and 1"). R&L don't adopt it; their
Eq. (11) lands `a ∈ [½, ⅔]` as an OUTCOME (⅔ at the innermost cell, `r_i = 0`).

## ⭐ Unread primary: Grant (1968)

I. P. Grant, "Numerical analysis of discrete ordinate methods," *J. Comput. Phys.*
**2**(4), 381–402 (1968), DOI `10.1016/0021-9991(68)90044-2` — ✅ CrossRef-verified.
**NOT in `scratch/literature/`.** It is (a) the origin of the weighted-diamond node
ansatz Eqs. (3)/(4) — the closure ORPHEUS ships, (b) the origin of the `[½,1]`
weight interval, (c) an alternative α determination ("assumes an explicit form for
α"). ⚠ OA: `is_oa False / closed`, no repo copy — a real acquisition decision.
