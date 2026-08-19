# DSA literature memo — consistent Diffusion Synthetic Acceleration for the ORPHEUS SN solver (issue #2, Phase 3)

**Purpose**: the literature record a future implementing session needs to build
consistent DSA for the slab DD SN solver and gate it against published spectral-radius
predictions. Every claim below carries source + equation number + printed page (and PDF
page for `Read pages=`). Extracted 2026-07-26 by literature-researcher; all equations
spot-verified against the rendered PDF pages (the scan is the SSOT).

---

## 0. Source inventory *(rev. 2, Adams-Martin/WWM — 2026-07-26)*

Seven primary sources are LOCAL in `scratch/literature/`. #1–#5 carry a Mistral-OCR
sidecar in `scratch/literature_ocr/<stem>.md` (`## p. N` = 1-based PDF page); #6–#7
arrived while the Mistral key was returning 401s (key worked at 00:50, dead by ~17:00
the same day) and carry a locally-generated pypdf **`.textlayer.md`** file instead
(same `## p. N` convention; equations verified on rendered pages). Regenerate their
true sidecars when the key is fixed: `tools/ocr_literature.py --glob 'Adams-Martin*'`
/ `--glob 'Warsa*'` (cache-idempotent).

| # | Paper | Journal | DOI | Local status | Page map (printed → PDF) |
|---|---|---|---|---|---|
| 1 | Alcouffe (1977), "Diffusion Synthetic Acceleration Methods for the Diamond-Differenced Discrete-Ordinates Equations" | NSE **64**(2), 344–355 | 10.13182/NSE77-1 | PDF + sidecar ✓ | printed ≈ PDF + 342 (printed 344 = PDF 2) |
| 2 | Larsen (1982), "Unconditionally Stable DSA Methods for the Slab Geometry Discrete Ordinates Equations. Part I: Theory" | NSE **82**(1), 47–63 | 10.13182/NSE82-1 | PDF + sidecar ✓ | printed ≈ PDF + 45 (printed 47 = PDF 2) |
| 3 | McCoy & Larsen (1982), "… Part II: Numerical Results" | NSE **82**(1), 64–70 | 10.13182/NSE82-2 | PDF + sidecar ✓ | printed ≈ PDF + 62 (printed 64 = PDF 2) |
| 4 | Adams & Larsen (2002), "Fast Iterative Methods for Discrete-Ordinates Particle Transport Calculations" | PNE **40**(1), 3–159 | 10.1016/S0149-1970(01)00023-3 | PDF + sidecar ✓ + **FULL prior extraction** `.claude/plans/archive/phase_i_survey_adams_larsen_2002.md` | journal p. N = PDF N − 2 |
| 5 | Morel (1982), "A Synthetic Acceleration Method for Discrete Ordinates Calculations with Highly Anisotropic Scattering" | NSE **82**(1), 34–46 | 10.13182/NSE82-A19026 | PDF + sidecar ✓ | printed ≈ PDF + 32 (printed 34 = PDF 2) |
| 6 | Adams & Martin (1992), "Diffusion Synthetic Acceleration of Discontinuous Finite Element Transport Iterations" — **the M4S paper** | NSE **111**(2), 145–167 | 10.13182/NSE92-A23930 | PDF ✓ + `.textlayer.md` (no Mistral sidecar yet) — **extracted, §7** | printed ≈ PDF + 143 (printed 145 = PDF 2) |
| 7 | Warsa, Wareing, Morel, McGhee & Lehoucq (2004), "Krylov Subspace Iterations for Deterministic k-Eigenvalue Calculations" — **the k-eigenvalue IRAM paper, NOT the degraded-DSA paper** | NSE **147**(1), 26–42 | 10.13182/NSE04-1 | PDF ✓ + `.textlayer.md` (no Mistral sidecar yet) — **extracted, §8** | printed ≈ PDF + 24 (printed 26 = PDF 2) |

Note that #2, #3, and #5 all sit in the same NSE issue (82(1), 1982) — the DSA issue.

**ABSENT from the local folder** (NOT acquired online per the delegation rule — user
to add if wanted):

- **Warsa, Wareing & Morel (2004)**, "Krylov Iterative Methods and the Degraded
  Effectiveness of Diffusion Synthetic Acceleration for Multidimensional S_N
  Calculations in Problems with Material Discontinuities," NSE **147**(2), 218–248 —
  the degraded-DSA robustness paper ("fully consistent DSA can degrade in multi-D
  heterogeneous problems, Krylov fixes it"). ⚠ Same volume/year/lead-authors as
  local #7 — do not conflate (see §8's identity block). Local #7's Ref. [25] cites it
  in submitted form under the early title "Diffusion Synthetic Acceleration—Part I:
  Deficiencies in Multi-Dimensional Heterogeneous Problems". A&L 2002 carries only
  the precursor thread (consistent-DSA heterogeneity degradation to σ ≈ 0.88 with
  "no known remedy", pp. 69/139, and the 2001 M&C proceedings version [300]).
- **Yavuz & Larsen (1988)**, Trans. ANS **56**, 305 — the reflective-BC DSA analysis
  (A&L [135], bibliography-only there; see §4.3/§6.6.1).

Bonus local holding (not part of this brief, relevant to the wider acceleration
program): `Hammer-Morel-Wang(2019)` nonlinear diffusion acceleration of least-squares
transport in void geometries.

---

## 1. Alcouffe 1977 — the founding consistent-DSA paper (DD SN)

R.E. Alcouffe, NSE **64**(2), 344–355 (1977). All printed pages p.NNN below; PDF page =
printed − 342. Equations (1)–(6) verified on PDF pp. 3–4; (12)–(23b) on PDF pp. 6–7;
(24)–(34) on PDF p. 8; Tables I–IV on PDF pp. 12–13.

### 1.1 The three continuous variants (§II.A, pp. 345–347)

One inner-iteration cycle: sweep the transport equation with the previous accelerated
scalar flux in the scattering source, then solve a *corrected diffusion equation* for the
next scalar flux. Tilde marks transport-iterate quantities (the modern half-iterate
φ^{ℓ+1/2}); untilded φ is the accelerated (diffusion-corrected) flux.

- **Sweep** (Eq. 1, p. 345): Ω·∇ψ̃_g^l + σ_g ψ̃_g^l = σ_{s,g→g} φ_g^{l−1} + QQ_g
  (isotropic scattering shown; the paper states the schemes extend to anisotropic).
- **Source-correction scheme** (linear; fixed-source problems). Corrected diffusion
  (Eq. 2, p. 345): −∇·D_g∇φ_g^l + σ_{Rg} φ_g^l = QQ_g **− R_g^l**, with
  D_g = 1/(3σ_tr), σ_{Rg} = σ_g − σ_{s,g→g}, and (Eq. 3, p. 346)

      R_g^l = ∇·J̃_g^l + ∇·D_g∇φ̃_g^l ,     J̃ = ∫ Ω ψ̃ dΩ, φ̃ = ∫ ψ̃ dΩ (Eqs. 4a/4b).

  R is the *diffusion defect of the transport iterate*: at convergence (φ = φ̃) the
  ∇·D∇ terms cancel and the transport balance ∇·J̃ + σ_R φ̃ = QQ remains (p. 346) —
  the accelerated answer IS the SN answer. R → 0 wherever the iterate satisfies Fick's
  law, so the low-order solve is exact precisely where SI is slow (diffusive regions).
- **Outer acceleration** (Eq. 5, p. 346): a multigroup corrected-diffusion solve with
  fission (χ_g Σ_{g′} νσ_{f,g′} φ_{g′}^{k+1}) and full group-transfer (Σ_{g′≠g}
  σ_{s,g′→g} φ_{g′}^{k+1}) held IMPLICIT at the new outer index, source QQ_g − R_g^k with
  the group-converged R. Strategy (p. 346): one inner per group per outer until the
  multigroup source converges, then converge the group fluxes; "the diffusion equation
  performs the outer iteration … the transport equation [supplies] the necessary
  corrections."
- **Diffusion-coefficient correction** (nonlinear; eigenvalue problems). Making Eq. (5)
  homogeneous by moving the correction INTO D (Eq. 6, p. 346):
  {D_g}_ii = −J̃_{gi}/∇_iφ̃_g (per-axis diagonal tensor) ⟹ R ≡ 0; inner (Eq. 7) and
  k_eff outer (Eq. 8, p. 347) use this D̃. This is the direct ancestor of the
  CMFD/NDA-family closures. Hazard (p. 347): D̃ can come out infinite or negative.
- **Removal correction** (nonlinear; the fallback): fold R into the removal instead
  (Eq. 9, p. 347): σ̃_{Rg}^k = σ_{Rg} + R_g^k/φ̃_g^k; inner (Eq. 10), outer (Eq. 11).
  Production practice (§II.B.5, p. 351, ONETRAN-DA): source correction for fixed-source;
  D-correction for eigenvalue; switch to removal correction cell-by-cell where D̃ < 0 or
  huge; if the removal also goes badly negative, drop acceleration for that group
  (p. 354). Fixed-source vs eigenvalue split: the source scheme NEEDS an inhomogeneous
  source (subcritical); the two nonlinear schemes exist to make the low-order problem
  homogeneous for k_eff.

### 1.2 The discrete DD SN equations (§II.B.1, p. 348)

Slab, one group (inner iteration suffices), cell centers i, edges i±1/2,
h_i = x_{i+1/2} − x_{i−1/2}:

- **DD SN** (Eq. 12): μ_m(ψ̃^l_{m,i+1/2} − ψ̃^l_{m,i−1/2}) + σ_i h_i ψ̃^l_{mi}
  = σ_{si} h_i φ^{l−1}_{0i} + QQ_i h_i.
- **Moments** (Eq. 13): φ_{ni} ≈ Σ_m w_m P_n(μ_m) ψ_{mi} with **Σ_m w_m = 1**
  (weights normalized to ONE — flag for the ORPHEUS quadrature map, §6.5).
- **Diamond closure** (Eq. 14): ψ̃^l_{mi} = ½(ψ̃^l_{m,i+1/2} + ψ̃^l_{m,i−1/2}).
- **P0 moment of (12)** (Eq. 15): φ̃^l_{1,i+1/2} − φ̃^l_{1,i−1/2} + σ_i h_i φ̃^l_{0i}
  = σ_{si} h_i φ^{l−1}_{0i} + QQ_i h_i.
- **P1 moment of (12)** (Eq. 16): **⅓**(φ̃^l_{0,i+1/2} − φ̃^l_{0,i−1/2})
  + **⅔**(φ̃^l_{2,i+1/2} − φ̃^l_{2,i−1/2}) + σ_i h_i φ̃^l_{1i} = 0.
  ⚠ OCR hazard confirmed live: the sidecar renders the ⅓ as ½; the printed page
  (PDF p. 6) clearly shows ⅓ (= Σw_mμ_m² for Σw=1, via μ² = ⅓P₀ + ⅔P₂). Trust ⅓.

Equation (16) is the structural heart: the DD P1-moment equation couples the
**cell-center current** φ_{1i} to **cell-edge scalar-flux differences** — a staggered
Fick's law. Every consistency statement below flows from honoring that staggering.

### 1.3 The inconsistent (Reed) scheme and why it fails (§II.B.2, pp. 348–349)

The "heretofore used" differencing (Reed, NSE **45**, 245 — Alcouffe's footnote 5 prints
year 1973; the A&L 2002 bibliography [51] and NSE volume numbering say **1971**): a
conventional **cell-centered** diffusion discretization. R discretized by integrating
Eq. (3) over cell i (Eq. 17, p. 348): edge currents φ̃_{1,i±1/2} plus a three-point
∇·D∇φ̃ built on CELL-CENTER fluxes φ̃_{0,i−1}, φ̃_{0,i}, φ̃_{0,i+1} with interpolated
edge coefficients D_{i+1/2} = ½(h_iD_i + h_{i+1}D_{i+1})/h_{i+1}; the corrected
diffusion equation (Eq. 18) is the matching cell-centered three-point scheme with
source QQ_i h_i − R_i^l. Compatible at convergence with the transport balance (Eq. 19):
φ̃_{1,i+1/2} − φ̃_{1,i−1/2} + σ_{Ri} h_i φ̃_{0i} = QQ_i h_i.

**The failure** (p. 348, verified on the page): the divergent problems have optically
thick regions, c = σ_s/σ near 1, mesh > 1 mfp. Diagnosis: in a homogeneous diffusive
region R_i of Eq. (17) does NOT → 0 unless h_i → 0, because "setting the second moment
terms in Eq. (16) to zero (the diffusion limit), a relationship for the first moment at
the cell centers results, while Eq. (17) utilizes the first moment at the cell
boundaries" (p. 348). The discrete diffusion operator disagrees with the discrete-SN
diffusion limit by O(1) at finite h ⟹ the correction source stays finite where it should
vanish ⟹ instability as Σ_t h grows.

### 1.4 The consistent scheme (§II.B.2, p. 349)

Derivation route (the paper's own, all on p. 349): (i) add the cell-i and cell-(i+1)
balance equations and use the P1 moment of the diamond relation to get the
**edge-centered balance** (Eq. 20):

    φ̃_{1,i+1} − φ̃_{1,i} + ½(σ_{R,i+1}h_{i+1}φ̃_{0,i+1} + σ_{Ri}h_iφ̃_{0i})
      = ½(QQ_i h_i + QQ_{i+1}h_{i+1}) ;

(ii) solve the converged Eq. (16) for the cell-center currents φ_{1i}, φ_{1,i+1} and
substitute ⟹ the exact edge-centered "diffusion form" of DD SN (Eq. 21), a three-point
relation in the EDGE fluxes φ_{0,i−1/2}, φ_{0,i+1/2}, φ_{0,i+3/2} with the exact
generalized coefficient

    D′_i = (1/3σ_i)·(1 + 2(φ_{2,i+1/2} − φ_{2,i−1/2})/(φ_{0,i+1/2} − φ_{0,i−1/2})) ;

when φ_2 ≡ 0, Eq. (21) IS the conventional edge-differenced diffusion equation.
(iii) Therefore difference R on the staggered (edge-centered) cell (Eq. 22):

    R_i^l = φ̃^l_{1,i+1} − φ̃^l_{1,i}
            + D_{i+1}(φ̃^l_{0,i+3/2} − φ̃^l_{0,i+1/2})/h_{i+1}
            − D_i(φ̃^l_{0,i+1/2} − φ̃^l_{0,i−1/2})/h_i ,     D_i = 1/(3σ_i),

— cell-CENTER currents, cell-EDGE scalar fluxes — and (iv) solve the edge-centered
corrected diffusion equation (Eq. 23, unknowns = edge fluxes φ^l_{0,i+1/2},
i = 1 … I−1 interior edges):

    −D_{i+1}(φ^l_{0,i+3/2} − φ^l_{0,i+1/2})/h_{i+1}
    + D_i(φ^l_{0,i+1/2} − φ^l_{0,i−1/2})/h_i
    + (σ̃_R h)^l_{i+1/2} φ^l_{0,i+1/2} = (QQh)_{i+1/2} − R_i^l ,

with the **nonlinear removal** (Eq. 23a):
(σ̃_R h)^l_{i+1/2} ≡ ½(σ_{R,i+1}h_{i+1}φ̃^l_{0,i+1} + σ_{Ri}h_iφ̃^l_{0i})/φ̃^l_{0,i+1/2}
and (QQh)_{i+1/2} = ½(QQ_i h_i + QQ_{i+1}h_{i+1}) (Eq. 23b). At convergence Eq. (20)
is recovered exactly, so the fixed point is the DD SN solution. The removal is a
**three-point (two-cell) weighted term**, not a one-point σ_R·φ_edge — this and the
edge-centered stencil are exactly the "differences from the conventional
discretization" A&L 2002 later names (their p. 56).

**Why nonlinear**: "to achieve this objective without assuming that the
diamond-difference relationship holds, we have introduced a nonlinearity into the
iteration scheme [Eq. (23a)]" (p. 349) — hence the abstract's "linear in its continuous
form and nonlinear in its differenced form." If one DOES assume diamond
(φ̃_{0i} = ½(φ̃_{0,i+1/2} + φ̃_{0,i−1/2})), (23a) linearizes to fixed ¼-weights —
that is the form analyzed in §II.B.3 and the form algebraically equivalent to Larsen's
four-step DD instance (A&L 2002, p. 54).

### 1.5 ⚠ Sign errata (load-bearing for any implementation from this paper)

The printed discrete pairs carry two typesetting slips, resolved unambiguously by
(a) the continuous scheme (2)–(3) ["QQ **−** R", R = +(∇·J̃ + ∇·D∇φ̃)], (b) each pair's
stated convergence identity [(17)/(18) ⟹ (19); (22)/(23) ⟹ (20)], and (c) the
model-problem restatement (24)–(26), which is fully sign-consistent:

| Printed | Slip | Correct reading |
|---|---|---|
| Eq. (17), p. 348 | leading signs "−φ̃_{1,i+1/2} + φ̃_{1,i−1/2} − D_{i+1/2}(…) + D_{i−1/2}(…)" define R as the NEGATIVE of discretized Eq. (3), yet Eq. (18) subtracts it | R_i = +[φ̃_{1,i+1/2} − φ̃_{1,i−1/2}] + [D_{i+1/2}(…) − D_{i−1/2}(…)] per Eq. (3); Eq. (18)'s "QQ_i h_i − R_i^l" then reduces to (19) at convergence |
| Eq. (23), p. 349 | prints "= (QQh)_{i+1/2} **+** R_i^l" (and drops the prime on R′) | "= (QQh)_{i+1/2} **−** R_i^l" with R per Eq. (22); only then does convergence give (20). Confirmed by Eq. (25) which prints "QQ_i − R_i^l" for the identical scheme |

Rule for the implementer: fix R by Eq. (3)/(22)/(26) (defect = ∇·J̃ + ∇·D∇φ̃, positive
sign) and SUBTRACT it from the low-order source, so the discrete ∇·D∇ terms cancel at
convergence leaving the transport balance. Never copy the printed (17)/(23) signs.

### 1.6 Stability demonstration (§II.B.3, pp. 349–350)

Model problem: slab, isotropic scattering, uniform h, constant σ, infinite medium (for
the N-matrix), diamond assumed in (23a) ⟹ linear scheme (24)–(26). Uniform-h forms
(all divided by h): sweep (24); corrected diffusion (25) with removal
¼σ_R(φ_{0,i+3/2} + 2φ_{0,i+1/2} + φ_{0,i−1/2}) and source QQ_i − R_i^l; and (26)
R_i^l = (φ̃^l_{1,i+1} − φ̃^l_{1,i})/h + (D/h²)(φ̃^l_{0,i+3/2} − 2φ̃^l_{0,i+1/2} +
φ̃^l_{0,i−1/2}).

The pivotal identity (Eq. 27, p. 350): taking the P1 moment of the sweep and inserting
it into (26),

    R_i^l = −(2D/h²)(φ̃^l_{2,i+3/2} − 2φ̃^l_{2,i+1/2} + φ̃^l_{2,i−1/2})

— the correction source is EXACTLY a second difference of the second Legendre moment of
the transport iterate. R ≡ 0 whenever φ̃_2 ≡ 0: the discrete scheme inherits the
continuous property "the low-order equation is exact in the P1 limit" at ANY h. This is
the cleanest single-equation statement of discrete consistency in the paper.

Error-iteration algebra (source set to 0): [**D** + (1−c)σ]φ̃_0^l = −2**D**φ̃_2^l
(Eq. 28, **D** = the edge diffusion matrix); φ̃_2^l = −**N**φ_0^{l−1} (Eq. 29, infinite
medium); iteration matrix **M**_s = 2[**D** + (1−c)σ]⁻¹**DN** (Eq. 30). With Reed's
eigenfunctions φ_k = cos((i+½)π/k): eigenvalues of **N** (Eq. 31)
λ_k = c⟨r_m²P₂(μ_m)/(r_m² + τ_k)⟩, r_m = σh/(2|μ_m|), ⟨f⟩ = Σw_mf_m,
τ_k = (1−cos π/k)/(1+cos π/k) (Eq. 32); of **D** (Eq. 33) η_k = (2D/h²)(1−cos π/k);
combined θ_k = 2η_kλ_k/[η_k + ½(1−c)σ(1+cos π/k)]. Bounds: λ_k < c/2 (using
⟨μ_m²⟩ = ⅓), and (Eq. 34, p. 350)

    θ_k < c / [1 + (1−c)(σh²/4D)·cot²(π/2k)] < c   for 0 ≤ c ≤ 1,

a STRICT inequality independent of σh²/4D ⟹ **spectral radius < c ≤ 1 for every mesh
size: unconditionally convergent** for this model problem. (Note the bound is ρ < c —
Alcouffe proves unconditional *stability*; the sharper ρ ≈ 0.2247c efficiency numbers
come from Larsen's Fourier analysis, §2 below. Note also λ_k is bounded by the FLATTEST
mode k → ∞, τ_∞ = 0.)

Contrast stated for Reed's scheme: with the cell-centered differencing, convergence
"depends upon the spatial mesh spacing" (p. 349, referring to Reed's demonstration).

### 1.7 Boundary conditions — what the paper does and does not give

- §II (p. 345): "We do not dwell on boundary conditions used with the transport and
  diffusion equations since they are the standard vacuum, reflective, and periodic
  conditions." The slab derivation presents INTERIOR rows only — Eq. (23) runs over
  interior edges i = 1 … I−1; the boundary-edge rows are never written.
- The implementation details (including curvilinear geometry and, implicitly, boundary
  rows) are deferred to Ref. 16, the ONETRAN-DA report (Alcouffe & Bosler, LASL, "to be
  published" — NOT in the local library).
- Empirical coverage: Table III's 2-D stability demonstration uses REFLECTIVE
  boundaries on two sides (left/bottom) with vacuum on the others, stable to 15 mfp
  cells — evidence, not derivation.
- ⟹ For the exact consistent boundary rows, use Larsen Part I's four-step boundary
  treatment (§2.4 below), not Alcouffe.

### 1.8 Two-dimensional extension (§II.C, pp. 351–353)

x-y DD SN (Eq. 35) + double diamond (Eq. 36); spherical-harmonic moments Y_p^q
(normalization block, p. 351–352); three discrete moments (Y_0^0, Y_1^0, Y_1^1) of
(35) give the 2-D moment equations (Eqs. 37–39); solving (38)/(39) for the two
cell-center current components yields Fick-plus-second-moment forms (Eq. 40) with
F, G → 0 in the diffusion limit; adding the four balance equations of cells (i,j),
(i+1,j), (i,j+1), (i+1,j+1) + diamond gives the **corner-centered** balance (Eq. 41);
substitution produces the corner-based 5-point corrected diffusion equation (Eq. 42)
with corner-averaged removal (42b), source (42c), and correction R_ij (Eq. 43) built
from corner-edge fluxes and cell-center currents. "To generate the correction term,
the flux moments must be evaluated on the spatial mesh cell corners" (p. 353) —
TWOTRAN-II was altered to compute corner angular fluxes (x-y, r-z, r-θ). This is the
vertex-centered ancestor of the multi-D 9-point consistent stencil (A&L 2002 §IV.D.4).
An O(h²) lumping (unnumbered display after Eq. 43, p. 353) collapses the
2×2-cell edge terms onto corner fluxes.

### 1.9 Numerical evidence (§III, pp. 353–355)

Convergence criterion 10⁻⁴ pointwise on the scalar flux; cylindrical geometry unless
stated. DSA converged on ALL problems; abstract's headline: 2–10× better than
rebalance/Chebyshev, most pronounced in large scattering regions.

**Table I (1-D total iterations; a = divergence):**

| Problem | NE | CY | CMR | FMR | DSA | time ratio (DSA/best other) |
|---|---|---|---|---|---|---|
| 1. LMFBR (28g, S4) | 1114 | 2382 | 3133 | a | **154** | 0.18 |
| 2. HTGR (9g, upscatter, S4) | — | 12000 | 6400 | — | **124** | 0.08 |
| 3. TREAT 1 (20g, S4) | — | 5047 | a | a | **188** | 0.06 |
| 4. TREAT 2 (44g, 20 upscatter) | — | — | — | — | **663** | — |
| 5. CTR (46g, S6, P3) | — | 1204 | — | — | **329** | 0.26 |
| 6. SPHERE (42g, S8, P3, adjoint source) | 1879 | — | a | — | **456** | 0.33 |

Note rebalance DIVERGES on 4 of the 6 entries where tried — the practical motivation.
Problem 6 exercised the abandon-acceleration fallback (removal term too negative).

**Table II** (TREAT 1, k_eff by outer): diffusion 1.2346 (−0.83%) → outer 1: 1.2451
(+0.01%) → outer 2: 1.2449 (−0.005%) → outer 3: 1.2450 (0.0%). One sweep + corrected
diffusion ≈ transport-quality k at ~2× diffusion cost (the "improved diffusion" mode,
§II.B.5 p. 351).

**Table III (2-D, homogeneous square, PURE scatterer, 5×5 mesh, reflective
left/bottom, flat source, cell size varied via σ):**

| cell (mfp) | FMR | DSA |
|---|---|---|
| 1 | 12 | **5** |
| 5 | 112 | **6** |
| 15 | 293 | **6** |

Mesh-thickness-INDEPENDENT DSA counts (5→6→6) vs FMR's 12→112→293 blowup: the
consistency payoff in one table. This is a c = 1 problem (pure scatterer).

**Table IV (2-D eigenvalue, total iterations):** problems 1/2 = Table-III square with
1/2 groups; problem 3 = 9-group Pu-core cylindrical reactor (r = 7.62 cm core,
half-height 1.966 cm, 7.62 cm depleted-U reflector; XS in LA-6571-PR): FMR 48/143/149;
CMR —/—/369; DSA **17/24/80**. (These are the numbers A&L 2002 re-quote at their
p. 136.)

**The fixup warning** (§IV, p. 355): TWOTRAN's negative-flux fixup "violates the
diamond equations … the compatible diffusion and transport difference equations that we
labored so carefully to construct are no longer compatible" — the original statement of
the fixup-breaks-consistency hazard (A&L 2002 M2 item 10).

### 1.10 Notation map (Alcouffe → ORPHEUS)

- μ_m = the slab direction cosine = ORPHEUS `mu_x`; m = 1…M ordinates.
- **Σ_m w_m = 1** (p. 348). ⟹ his φ_0 = quadrature AVERAGE of ψ (½∫ψdμ in continuous
  terms), his ⟨μ²⟩ = ⅓. Check ORPHEUS's slab quadrature normalization (Σw = 2
  convention would rescale all moment equations by 2) before transcribing (13)–(16).
- σ = σ_t (total); σ_R = removal = σ_t − σ_{s,g→g}; c = σ_s/σ; D = 1/(3σ_tr)
  (transport-corrected σ possible; = 1/(3σ_t) for isotropic).
- Superscript l = inner iterate; k = outer iterate; tilde = post-sweep transport
  quantity (≡ the modern half-iterate ℓ+1/2); untilded = accelerated.
- His low-order solves for the FULL new flux with defect-corrected source
  (QQ − R), NOT for an additive correction δφ. Equivalent to the additive
  correction form (Larsen's f-form, §2) when the removal is linear: subtract the
  edge-balance of the transport iterate from Eq. (23). ORPHEUS should implement ONE
  form and note the equivalence.
- Cell-edge unknowns: the consistent low-order system lives on EDGES (i+1/2), size
  I−1 interior + 2 boundary rows; the inconsistent one on centers.

## 2. Larsen 1982 Part I — the four-step recipe

E.W. Larsen, NSE **82**(1), 47–63 (1982). Printed p.NN; PDF = printed − 45. Equations
(10)–(28) verified on PDF pp. 6–9; (29)–(40) on PDF pp. 10–12; Table I on PDF p. 6.
Scope: slab SN with **weighted diamond (WD, includes DD), linear characteristic (LC),
linear discontinuous (LD), linear moments (LM)** schemes; one group; up to linearly
anisotropic scattering (with the n ≥ 2 generalization discussed); Gauss-Legendre
quadrature, **Σ_m ω_m = 1**. Iterate convention (p. 51): superscript l+1/2 = computed by
transport; l, l+1 = computed by diffusion — matches A&L 2002 and the ORPHEUS plan.

### 2.1 The continuous prototype (§I, pp. 48–50) — the essential idea

Unaccelerated: sweep (1) μ∂ψ^{l+1/2}/∂x + ψ^{l+1/2} = cφ₀^l + S (mfp units, σ_T = 1);
update (2) φ₀^{l+1} = φ₀^{l+1/2} = ½∫ψ^{l+1/2}dμ. Fourier ⟹ ω = (c/2)∫dμ/(1+λ²μ²),
spr = c; the slow modes are λ ≈ 0, "nearly linear functions of μ" (p. 48).

Moments of the sweep (3a/3b), φ_n = ½∫P_nψdμ:

    dφ₁^{l+1/2}/dx + φ₀^{l+1/2} = cφ₀^l + S ,
    (2/3)dφ₂^{l+1/2}/dx + (1/3)dφ₀^{l+1/2}/dx + φ₁^{l+1/2} = 0 .

**The essential idea** (p. 49, the paper's own words: "the derivation of Eqs. (5) from
Eqs. (3) is the essential idea in this paper"): promote the iteration index to l+1 on
every φ₀, φ₁ slot; lag φ₂ at l+1/2 (Eqs. 5a/5b):

    dφ₁^{l+1}/dx + (1−c)φ₀^{l+1} = S ,
    (2/3)dφ₂^{l+1/2}/dx + (1/3)dφ₀^{l+1}/dx + φ₁^{l+1} = 0 .

If ψ is linear in μ (⟹ φ₂ = 0), the promoted system determines φ₀^{l+1}, φ₁^{l+1}
EXACTLY in one solve — the accelerator is exact precisely on the slow (λ ≈ 0,
linear-in-μ) modes. Subtracting (3) from (5) gives the equivalent **correction (f)
form** (6a/6b):

    −(1/3)d²f₀^{l+1}/dx² + (1−c)f₀^{l+1} = c(φ₀^{l+1/2} − φ₀^l) ,
    φ₀^{l+1} = φ₀^{l+1/2} + f₀^{l+1} ,

with f_n^{l+1} = φ_n^{l+1} − φ_n^{l+1/2}. Gelbard-Hageman result (Eq. 7, p. 49):
spr = sup_λ [cλ²/(λ²+3(1−c))]·|∫₋₁¹ P₂(μ)/(1+λ²μ²)dμ| < **0.23c** (bound computed
numerically). Alcouffe's Eq. (8) full-flux form (with the −d/dx[⅓dφ₀^{l+1/2}/dx +
φ₁^{l+1/2}] source) is "algebraically equivalent" (p. 50).

**History as Larsen tells it** (pp. 49–50): Reed's discretized synthetic method
(DD transport + central-differenced diffusion) has spr ≥ c/[4/(3h²) + 1 − c] ⟹
**unstable (spr > 1) for h > h\* = 2/[3(2c−1)]^{1/2}**; c ≈ 1 ⟹ h\* ≈ 2/√3 ≈ 1.15 mfp.
McCoy showed Gelbard-Hageman's 5-point vertex scheme reduces in slab to an edge scheme
with EXACTLY Reed's spectral radius (p. 50) — so both classic discretizations die at
thick cells. Alcouffe fixed it for DD; refined analyses [16,17] sharpened his spr < c
to **spr ≤ 0.23c for all h**. And the transfer warning (p. 50): the DD-consistent
diffusion discretization "does not give unconditionally stable acceleration for the
weighted diamond (WD) differenced transport equation" [17,24] — consistency is
per-scheme, not a fixed low-order recipe.

### 2.2 The four-step procedure, discrete (§II, pp. 51–54) — WD instance

The recipe ORPHEUS's 3a computational derivation must reproduce. Steps, with the WD
equations:

**Step 1 — the discretized transport iteration.** Balance (10a) + WD closure (10b):

    (μ_m/h_i)(ψ^{l+1/2}_{m,i+1/2} − ψ^{l+1/2}_{m,i−1/2}) + σ_{Ti}ψ^{l+1/2}_{mi}
        = σ_{S0i}φ₀ᵢ^l + 3σ_{S1i}μ_mφ₁ᵢ^l + S_{mi} ,
    ψ^{l+1/2}_{mi} = ((1+α_{mi})/2)ψ^{l+1/2}_{m,i+1/2} + ((1−α_{mi})/2)ψ^{l+1/2}_{m,i−1/2} ,

weights: |α_{mi}| ≤ 1, μ_mα_{mi} ≥ 0, antisymmetric α(−μ_m) = −α(μ_m) (12a–c).
Members: DD α = 0; fixed-weight WD α = a_i·sgn(μ_m); step α = sgn(μ_m); step
characteristic α = coth(ε_{mi}) − 1/ε_{mi}, ε_{mi} = σ_{Ti}h_i/(2μ_m).

**Step 2 — take the L₀ and L₁ moments of BOTH equations** (L_n· = Σ_m P_n(μ_m)·ω_m).
Definitions (14a–c): ρ_i = L₁α_i = Σμ_mα_{mi}ω_m; γ_{mi} = α_{mi} − 3ρ_iμ_m;
β_{mi} = μ_mα_{mi} − ρ_i (so L₀γ = L₁γ = L₂γ = 0, L₀β = L₁β = 0; for linear-in-μ ψ:
L₀γψ = L₀βψ = φ₂ = 0, Eq. 15). Result = FOUR moment equations (16a–d): the two balance
moments (16a: P0; 16b: P1 with the (2/3h)Δφ₂ + (1/3h)Δφ₀ + σ_Tφ₁ structure) AND the
two closure moments (16c: φ₀ᵢ = ½(edge sum) + (3/2)ρ_i(φ₁ edge diff) + ½L₀γ_i(ψ edge
diff); 16d: φ₁ᵢ = ½(edge sum) + ½ρ_i(φ₀ edge diff) + ½L₀β_i(ψ edge diff)). Taking
moments of the CLOSURE, not just the balance, is what distinguishes the four-step from
naive P1-of-the-balance — the closure moments carry the scheme's α-weights into the
low-order operator via ρ_i.

**Step 3 — promote the iteration index** (the "discrete P1 approximation" slot
assignment): in (16a–d) replace l+1/2 → l+1 on every φ₀, φ₁ (edge and cell-average) and
l → l+1 in the scattering terms; LAG at l+1/2: φ₂ and the residual angular terms L₀γψ,
L₀βψ (Eqs. 18a–d, p. 53). If ψ^{l+1/2} were linear in μ, (18) determines the new
moments exactly (p. 53).

**Step 4 — subtract & reduce to one tridiagonal.** Subtract (16) from (18); define the
EDGE corrections f^{l+1}_{n,i+1/2} = φ^{l+1}_{n,i+1/2} − φ^{l+1/2}_{n,i+1/2}, n = 0,1
(19). The φ₂ and L₀γψ/L₀βψ terms CANCEL in the subtraction ⟹ the correction system
(20)–(21) is closed in f₀, f₁ alone. Eliminate cell-average corrections (21) → the
two-equation discrete P1 system (22a/22b) → compact form (24a/24b) with coefficients
(23a–f):

    σ̂_{Ri} = (σ_{Ti} − σ_{S0i}) / [1 + (3/2)ρ_i(σ_{Ti} − σ_{S0i})h_i] ,
    σ̂_{Si} = σ_{S0i} / [1 + (3/2)ρ_i(σ_{Ti} − σ_{S0i})h_i] ,
    g₀ᵢ^{l+1/2} = σ̂_{Si}h_i(φ₀ᵢ^{l+1/2} − φ₀ᵢ^l)      [the DSA residual source] ,
    D_i = 1/[3(σ_{Ti} − σ_{S1i})] + (1/2)ρ_i h_i ,
    a_i = σ_{S1i}/(σ_{Ti} − σ_{S1i}) ,   g₁ᵢ^{l+1/2} = a_i(φ₁ᵢ^{l+1/2} − φ₁ᵢ^l) .

Write f₁ at edge i+1/2 from cell i (25) and from cell i+1 (26); eliminating f₁ gives
**the consistency condition = the discrete diffusion equation** (27), a tridiagonal
system in the edge corrections f₀:

    −(D_{i+1}/h_{i+1})(f₀_{i+3/2} − f₀_{i+1/2}) + (D_i/h_i)(f₀_{i+1/2} − f₀_{i−1/2})
    + ¼[σ̂_{R,i+1}h_{i+1}(f₀_{i+3/2} + f₀_{i+1/2}) + σ̂_{Ri}h_i(f₀_{i+1/2} + f₀_{i−1/2})]
    = ½(g₀_{i+1} + g₀ᵢ) − (g₁_{i+1} − g₁ᵢ) ,    i = 1 … I−1 ,

then update the cell averages (28a/28b):

    φ₀ᵢ^{l+1} = φ₀ᵢ^{l+1/2} + (½ − ¾ρ_iσ̂_{Ri}h_i)(f₀_{i+1/2} + f₀_{i−1/2}) + (3/2)ρ_i g₀ᵢ ,
    φ₁ᵢ^{l+1} = φ₁ᵢ^{l+1/2} + (½ρ_i − D_i/h_i)(f₀_{i+1/2} − f₀_{i−1/2}) + g₁ᵢ .

One full iteration = sweep (10) → tridiagonal (27) → updates (28). Structure notes:

- **DD instance** (α = 0 ⟹ ρ = 0): σ̂_R = σ_R, σ̂_S = σ_S0, D = 1/[3(σ_T−σ_S1)] — the
  operator reduces to Alcouffe's edge-centered scheme with the ¼ three-point removal;
  (27)'s LHS is then identical in form to Alcouffe Eq. (23)/(25). For ρ ≠ 0 the
  four-step generates *h-dependent* coefficients: D grows with ρh and σ̂_R is h-damped
  — consistency is NOT "use the analytic D and σ_R," it bakes the closure weights in.
- **Both accelerated moments**: φ₀ AND φ₁ are updated (Morel's idea imported, p. 50–51)
  — needed under anisotropic scattering (§5 below).
- Higher moments stay unaccelerated: φ_n^{l+1} = φ_n^{l+1/2}, n ≥ 2 (p. 54) — linear,
  adequate while σ_{Sn} ≪ σ_T for n ≥ 2; Morel's nonlinear alternative otherwise.

### 2.3 The f-form vs the full-flux form (p. 53–54) — implementation guidance

Larsen derives both and RECOMMENDS the f-form (10)+(27)+(28) over the "original"
(10)+(18):

1. **Storage**: 5 arrays (f₀ edge; g₀, g₁, φ₀^{l+1/2}, φ₁^{l+1/2}) vs 9.
2. **Fixup robustness**: at convergence (27) → "0 = 0" identically, whereas (18)
   reduces to the transport equations only if the closure (10b) holds in every cell —
   "If we use negative flux fixups … Eq. (10b) will not hold in cells where fixups are
   used, and for such cells, Eqs. (18) must be modified … the recommended equation,
   Eq. (27), does not need to be modified" (p. 54; credit to Aull). This resolves
   Alcouffe's §IV fixup-incompatibility worry STRUCTURALLY: the correction form is
   insensitive to where the transport iterate came from.

### 2.4 Boundary conditions (§III, pp. 56–57) — the exact rows

Two distinct BC families; the guiding principle (p. 56): "If the exact numerical
transport solution ψ_m is a linear function of μ_m, then the numerical process should
generate this solution in the first iteration."

**(a) The correction equation (27)** — since f₀ → 0 at convergence, "any homogeneous
boundary conditions (i.e., ones that admit the zero solution) are permissible" (p. 57);
but to preserve transport features, define the diffusion edge angular flux
Ψ^{l+1}_{m,edge} = φ₀^{l+1}_{edge} + 3μ_mφ₁^{l+1}_{edge} and impose the transport BC on
it. Resulting rows:

- Prescribed incident flux (vacuum = zero incident), left/right (38a/38b):
  0 = (Σ_{μ_m>0} μ_mω_m)·f₀,₁/₂ + ½f₁,₁/₂  [right: sum over μ_m < 0]. This is the
  DISCRETE Marshak condition: the half-range coefficient Σ_{μ>0}μ_mω_m (= A&L's γ_N,
  → ¼ as N → ∞ under Σω = 1) replaces the continuum ¼; using the continuum value at
  small N is a quiet consistency break (A&L 2002 M2 item 4).
- Reflecting (39a/39b): f₁,edge = 0.
- Periodic (40): f_{n,1/2} = f_{n,I+1/2}, n = 0, 1.
- Close the rows by eliminating the boundary-edge f₁ via the one-sided cell relations
  (25) at i = I and (26) at i = 0 (explicitly written for LC as (62a/62b), p. 59;
  the WD forms are the same with the WD coefficients).

**(b) The starting iterate** (l = 0 with ψ^{1/2} = 0): solve the CONSISTENT full-flux
diffusion system (30)–(32) — same operator as (27), source (Q₀_{i+1} + Q₀ᵢ) −
(Q₁_{i+1} − Q₁ᵢ) with Q₀ᵢ = h_i(L₀S_i)/[2 + 3ρ_i(σ_{Ti}−σ_{S0i})h_i], Q₁ᵢ =
L₁S_i/(σ_{Ti}−σ_{S1i}); updates (31a/31b). Its BCs (33)–(37): the linear-in-μ edge
ansatz ψ_{m,1/2} = φ₀,₁/₂ + 3μ_mφ₁,₁/₂ (33); incident-flux row (34a):
Σ_{μ>0}μ_mψ_{m,1/2}ω_m = (Σ_{μ>0}μ_mω_m)φ₀,₁/₂ + ½φ₁,₁/₂; reflecting (35): φ₁,edge = 0;
periodic (36); one-sided current closures (37a/37b). Note ½ = 3·Σ_{μ>0}μ_m²ω_m
(half-range second moment, exact for Gauss-Legendre with Σω = 1).

**S2 exactness** (p. 56): K_2 = 0 because with two ordinates every ψ_m IS linear in μ;
with the BCs (34)/(35)/(36) "these boundary conditions produce the exact transport
fluxes in one (diffusion) calculation" (p. 57) — S2-SN ≡ diffusion, discretely. A sharp
implementation test: an S2 run must converge in ONE accelerated iteration (to roundoff)
if and only if the BC rows are exactly consistent.

### 2.5 The Fourier analysis and the spectral-radius bounds (§III, pp. 55–56)

Model problem: infinite medium, isotropic (σ_S1 = 0), constant σ, uniform h, mfp units
(σ_T = 1, σ_S0 = c). Ansatz: φ₀ᵢ^l = exp(jλx_i); ψ^{l+1/2} edge/center amplitudes
t_m, u_m; f₀ amplitude v; φ₀ᵢ^{l+1} = ω·exp(jλx_i). Exact eigenvalue (Eq. 29, p. 55):

    ω = c·[η sin²ξ / (η sin²ξ + (1−c)cos²ξ)]
        ·[((2θ₁ + 3hθ₂)/(2 + 3hρ) + θ₃) / (1 + (3/2)ρh(1−c))] ,   ξ = λh/2 ,

    η = (4/3)[1/h + (3/2)ρ(1−c)][1/h + (3/2)ρ] ,
    θ₁ = Σ_m (1 − 3μ_m²)·cos²ξ/[cos²ξ + (α_m + 2μ_m/h)²sin²ξ]·ω_m ,
    θ₂ = Σ_m (ρ − μ_mα_m)·cos²ξ/[cos²ξ + (α_m + 2μ_m/h)²sin²ξ]·ω_m ,
    θ₃ = Σ_m (α_m − 3ρμ_m)·(α_m + 2μ_m/h)/[cos²ξ + (α_m + 2μ_m/h)²sin²ξ]·ω_m .

Since η > 0 and 0 ≤ c ≤ 1: |ω| ≤ c·|(2θ₁ + 3hθ₂)/(2 + 3hρ) + θ₃| (equality at c = 1)
⟹ **spr ≤ c·K_N**, K_N = sup over 0 < h < ∞, 0 ≤ ξ ≤ π/2, admissible α. Numerically
K_N is maximized at h = ∞, ξ = π/2, |α_m| = 1, giving the closed form (p. 56):

    K_N = 1 − 3ρ² ,   ρ = Σ_m |μ_m|ω_m  →  ½ as N → ∞  ⟹  K_∞ = ¼ .

**Table I (printed p. 51; the L1 gate values for ORPHEUS):**

| N | WD | LC | LD | LM |
|---|---|---|---|---|
| 2 | 0 | 0 | 0 | 0 |
| 4 | 0.185 | 0.236 | 0.238 | 0.223 |
| 6 | 0.220 | 0.281 | 0.279 | 0.260 |
| 8 | 0.233 | 0.296 | 0.291 | 0.271 |
| 12 | 0.243 | 0.308 | 0.298 | 0.278 |
| 16 | 0.246 | 0.312 | 0.299 | 0.280 |
| 20 | 0.247 | 0.314 | 0.300 | 0.282 |
| 24 | 0.248 | 0.315 | 0.300 | 0.282 |

All K_N < 1/3 ⟹ spr ≤ cK_N < c/3 for every scheme, every N, every h, every c ≤ 1:
**unconditional stability + effectiveness** (abstract's claim). Validity: model problem
only (infinite medium, isotropic, uniform); BUT since the method is linear and the model
has no boundaries, "convergence rates and stability for these problems cannot depend on
the choice of the initial iterate or boundary conditions" (p. 55). ⚠ These K_N are the
WD-family sup over |α_m| = 1 — for pure DD (α = 0) at FINITE h the spectral radius is
smaller; the DD-specific ρ(c, Σ_th) curve is in Part II (§3) and A&L 2002 Eq. (3.65).

### 2.6 LC / LD / LM arms (§IV–V, pp. 57–63) — for the future ORPHEUS LD arm

- **LC** (41)–(47): cell-wise analytic solution with linear scattering-source
  representation; closure = WD relation + slope term: ψ_{mi} = WD(edges) −
  α_m(σ_S/σ_T)ζᵢ^l, α_m = coth(ε_m) − 1/ε_m. Six unknowns (φ₀,φ₁,ζ₀,ζ₁,f₀,f₁); same
  four-step; final system (56)–(61): SAME tridiagonal shape, coefficients now
  D_i = 1/[3(σ_T−σ_S1)] + ½ρh(1 − σ_S0/σ_T), σ̂_R = σ_T(σ_T−σ_S0)/[σ_T +
  (3/2)ρh(σ_T−σ_S0)(σ_T−σ_S1)], sources g₀, g₁ picking up slope-difference terms
  (60)–(61); BC closures (62); starting iterate (63)–(66).
- **LD/LM** (67)–(88): jointly treated. LD = zeroth (70a) + first spatial balance (70b:
  (3μ_m/h)(ψ_{i+1/2} + ψ_{i−1/2} − 2ψ_i) + σ_Tξ_m = source-slope) + one-sided slope
  definition (68) ⟹ WD-form closure with **α_m = ε_m/(3 + |ε_m|)** (69a); LM uses the
  exact first-moment (46) ⟹ α_m = coth(ε_m) − 1/ε_m (same α as LC but different slope
  equation). Four-step ⟹ tridiagonal (78) + slope updates (79) + moment updates (80),
  coefficients (81a–d) (D_i and σ̂_R = (σ_T−σ_S0)/(1+Δ_i) with the (81d) Δ_i rational
  function), sources (82)–(83), BC closures (84), starting iterate (85)–(88).
- Meta-point for ORPHEUS architecture: in ALL four schemes the low-order operator is
  the SAME edge tridiagonal shape — only the rational-function coefficients
  (D_i, σ̂_R, σ̂_S, g's) change with the closure α. A computational moment-reduction
  (issue #2's 3a) should reproduce exactly this: operator shape fixed by the balance,
  coefficients generated by the closure.

### 2.7 Scope limits Part I itself states (§VI, p. 63)

- Slab only; "should be straightforwardly applicable to any differencing scheme in slab
  geometry possessing linear relationships between the cell-edge and cell-averaged
  quantities"; expected to extend to 1-D curvilinear (not done here).
- Multi-D: the four-step applies, but converting the moment equations to a single
  vertex-differenced diffusion equation is proven possible only for DD (Alcouffe's
  Ref. 11); other schemes "under investigation" (this is the gap the M4S/WWM line later
  fills — see §4).
- Eigenvalue problems: "currently examining" — Part I is fixed-source only.
- Flux-fixup WD weights (iteration-dependent α, violating (12c)) not covered by the
  analysis; generalization "under investigation."
- Nonlinear schemes (exponential) excluded — "not amenable to the analysis."

## 3. McCoy & Larsen 1982 Part II — measured stability envelopes

D.R. McCoy & E.W. Larsen, NSE **82**(1), 64–70 (1982). Printed p.NN; PDF = printed −
62. All tables transcribed from the rendered pages (PDF pp. 3–7); the sidecar carries
the prose only. Purpose: experimental confirmation that the Part-I model-problem bounds
(spr ≤ cK_N, K_2 = 0, K_N < 1/3 for N ≥ 4 — restated as Eqs. 1a–c, p. 64) survive
finite media, anisotropic scattering, variable σ and nonuniform mesh; plus the fixup
study. "K_N < 1/3 translates into … errors reduce by at worst one order of magnitude
for every two iterations, for any spatial mesh" (p. 64).

### 3.1 Test problems

- **Homogeneous slab**: 8 cm, 8 equal 1-cm cells, REFLECTING left / vacuum right,
  constant source in the left half, S4 Gauss-Legendre, pointwise 10⁻⁴, O(1) initial
  error; σ_T ∈ {1, 2, 4, 6, 10, 20} cm⁻¹ so cell optical width = σ_T mfp; c ∈ {0.98,
  0.8}. Unaccelerated DD needs > 40 (c = 0.8) and 200–300 (c = 0.98) iterations (p. 66).
- **Model shielding problem** (Fig. 1, p. 66; reflecting left / vacuum right; S4;
  linearly anisotropic; 10⁻⁴): four regions, cross sections in cm⁻¹:
  water 0–12 cm (σ_T = 3.3333, σ_S0 = 3.3136, σ_S1 = 0.9256, S = 1.0); water 12–15 cm
  (same, S = 0); iron 15–21 cm (σ_T = 1.3333, σ_S0 = 1.1077, σ_S1 = 0.0367, S = 0);
  water 21–30 cm (same as region 2). Five meshes, cells/region = (40,10,8,30),
  (20,5,4,15), (12,3,2,9), (4,1,2,3), (1,1,2,1) ⟹ per-region σ_Th = (1,1,1,1) up to
  (40,10,4,30) mfp. [Ready-made ORPHEUS regression problem: heterogeneous, anisotropic,
  reflective+vacuum, cells to 40 mfp.]
- **Ref. 5 model problem** (Alcouffe-Larsen-Miller-Wienke, NSE **71**, 111 (1979)):
  meshes M6/M3/M1 = 0.5/1.0/3.0 mfp cells, S8, 10⁻⁵ (used in the fixup study).

### 3.2 Consistent DSA without fixups — mesh-independence confirmed

**Table I (p. 65)** — WD transport + its OWN four-step DSA, c = 0.98 slab, iterations:

| σ_T (= cell mfp) | a=0.0 (DD) | a=0.1 | a=0.5 | a=1.0 |
|---|---|---|---|---|
| 1 | 4 | 4 | 4 | 4 |
| 2 | 4 | 4 | 3 | 4 |
| 4 | 5 | 5 | 4 | 4 |
| 6 | 7 | 5 | 4 | 4 |
| 10 | 5 | 6 | 5 | 4 |
| 20 | 4 | 5 | 6 | 5 |

"fully indicative of all our other results … a maximum of about eight iterations …
for any size of spatial mesh" (p. 65).

**Table III (p. 66)** — homogeneous slab, all four schemes with their own DSA:

| σ_T | c=0.98: DD / LC / LD / LM | c=0.8: DD / LC / LD / LM |
|---|---|---|
| 1 | 4 / 4 / 4 / 4 | 5 / 5 / 5 / 5 |
| 2 | 4 / 4 / 4 / 4 | 7 / 5 / 5 / 5 |
| 4 | 4 / 4 / 4 / 4 | 5 / 7 / 8 / 8 |
| 6 | 7 / 5 / 5 / 5 | 4 / 6 / 7 / 7 |
| 10 | 5 / 7 / 10 / 8 | 4 / 5 / 6 / 6 |
| 20 | 4 / 5 / 7 / 7 | 4 / 5 / 5 / 5 |

**Table IV (p. 66)** — shielding problem: DD 5–7, LC 5–7, LD 5–7, LM 5–6 across ALL
five meshes (up to 40-mfp cells). Contrast: **LD + coarse-mesh rebalance (ONETRAN)
needs 440 iterations** on this problem (p. 66).

### 3.3 The partial-consistency experiment — Table II (p. 65), the memo's sharpest result

"Diamond-synthetic" = WD transport (a ≠ 0) accelerated by the DD-derived diffusion
(Part-I Eqs. I(27), I(28), I(34)–(37) with ρ_i = 0 — i.e., the right OPERATOR SHAPE
with the wrong closure constant, including the BC rows):

| σ_T | a=0.0 | a=0.1 | a=0.5 | a=1.0 |
|---|---|---|---|---|
| 1 | 4 | 4 | 6 | 8 |
| 2 | 4 | 4 | 9 | 17 |
| 4 | 5 | 6 | 18 | 91 |
| 6 | 7 | 8 | 43 | div |
| 10 | 5 | 13 | div | div |
| 20 | 4 | 22 | div | div |

"for any nonzero value of a, no matter how small, we find experimentally that the
diamond-synthetic iteration scheme becomes unstable for sufficiently large meshes"
(p. 65). Two morals stated in §IV (pp. 69–70): (i) partial consistency has NO safe
dose — the degradation threshold just moves to thicker cells as the inconsistency
shrinks; (ii) **"it is essential that both the diffusion equation AND its boundary
conditions be compatible with the transport calculation"** — Table II's scheme was
incompatible in both, and the sensitivity is worst for highly anisotropic problems.
[This also models the fixup effect: a fixup makes the effective transport closure
depart from the one the diffusion operator was derived from, "which roughly
corresponds to the use of an increasing number of fixups," p. 65.]

### 3.4 P1 vs P0 acceleration under anisotropy — Table V (p. 67)

Modified shielding: water σ_S1 = βσ_S0, β ∈ {0.5, 0.8, 0.9, 0.99}; "P1" = accelerate
φ₀ AND φ₁ (Part I); "P0" = accelerate φ₀ only. Iterations (DD/LC/LD/LM):

| β | P1 fine | P1 coarse | P0 fine | P0 coarse |
|---|---|---|---|---|
| 0.5 | 5/5/5/5 | 5/6/5/6 | 9/9/9/9 | 7/9/9/7 |
| 0.8 | 5/5/5/5 | 5/6/6/6 | 19/18/18/18 | 8/20/17/12 |
| 0.9 | 5/5/5/5 | 5/6/7/6 | 29/28/28/28 | 8/31/23/16 |
| 0.99 | 5/5/5/5 | 5/6/6/7 | 64/64/65/64 | 35/73/47/22 |

P1 acceleration is β-INDEPENDENT (5–7 always); P0-only degrades roughly like the
Morel ρ = σ_S1/σ_T prediction (§5) — at β = 0.99, 64 vs 5 iterations. "for highly
anisotropic materials, it is essential to accelerate both the scalar fluxes and
currents" (p. 66); verified here for LC/LD/LM (Morel had shown DD).

### 3.5 Fixup study (§III, pp. 67–69) — Tables VI–IX

Fixup locations: (1) starting-calculation source can go negative "because of the
three-point removal term in the discretized diffusion equations" (p. 67) — set
negative cell-average source (and slope) to zero; (2) transport fixups — DD
set-to-zero; LC/LD/LM slope-limiting (minimum slope change making the linear
representation nonnegative; LD's fixup also alters the EXITING flux — "hard" fixup —
which LC/LM avoid); (3) diffusion output can go negative — if any cell-average
transport source < 0, DISCARD the diffusion result and use the previous fixed-up
transport source (acceleration off for that iteration).

Headline numbers (a = negative cell-average flux would occur without fixup;
b = divergence):

- **Table VI (c = 0.98 slab, p. 68)**: DD fine → 4–10, σ_T = 8: 40^a, σ_T ≥ 9: a,b
  (divergence). ONEDANT's nonlinear DD variant: same shape (39^a at 8, then a,b).
  LC/LD/LM: fine 4–8; a window of divergence (LC at 9; LD at 10; LM at 11) then
  RE-convergence at still-larger meshes but at unaccelerated-like counts
  (LC 269–452^a, LD 124–442^a, LM 166–358^a).
- **Table VII (c = 0.8, p. 68)**: DD dies from σ_T = 4 up; LC converges but degraded
  (36–68^a); LD diverges at σ_T ≥ 6; LM 38–59^a.
- **Table VIII (shielding, p. 68)**: DD (both variants) diverges on the three coarse
  meshes; LC/LD/LM converge on ALL meshes (LC 5–7, LM 5–7; LD 27 and 22 on two
  meshes — the hard LD fixup interfering via a negative STARTING solution, p. 68).
- **Table IX (Ref. 5 problem, S8, p. 69)**: M6/M3 all 7; M1 (3-mfp cells): DD a,b;
  LD a,b; LC 8; LM 8. [LD + rebalance in ONETRAN: 50 iterations on M6.]

Summary (p. 69): "the application of a negative flux fixup on a spatial mesh for which
a method without fixup would yield (or would almost yield) negative cell-averaged
fluxes is an anathema to all of the acceleration methods." Remedies (p. 70): refine
the mesh (source problems); turn the fixup OFF (eigenvalue problems); long-term —
positive-by-construction schemes or fixup-compatible acceleration ("difficult and
long-standing problem").

### 3.6 Cost model and multigroup status (§IV, p. 69)

- CPU: diffusion-solve/transport-sweep time ratio **< 0.1 at S4**, smaller at higher
  N ⟹ "the CPU time required to perform one accelerated iteration is essentially equal
  to … one unaccelerated iteration" — the 1-D cost model for ORPHEUS's gate economics
  (contrast A&L 2002's multi-D estimate: diffusion solve ≈ one S4 sweep).
- Multigroup: NOT tested; one-group results = inner-iteration acceleration; the outer
  extension is "conceptually straightforward" citing the DD experience
  (ONEDANT/Alcouffe).
- Morel's σ_S2 caveat (p. 69): P1 (φ₀+φ₁) acceleration itself degrades as σ_S2 → σ_S0,
  "very similar to the case of P0 acceleration [degrading] as σ_S1 approaches σ_S0";
  Morel's remedy = after the diffusion solve, rescale ALL EVEN transport angular-flux
  moments by one factor and ALL ODD moments by another so the new φ₀, φ₁ equal the
  diffusion values (nonlinear; see §5).

## 4. Adams & Larsen 2002 — the consistency taxonomy and everything around it

M.L. Adams & E.W. Larsen, PNE **40**(1), 3–159 (2002). Journal p. N = PDF page N−2.
A **complete** prior extraction exists: `.claude/plans/archive/phase_i_survey_adams_larsen_2002.md`
(equations transcribed from page images; OCR hazard ⅓Λ² → "½A2" documented there).
This section distills the DSA-relevant core; the survey holds the rest (Krylov history,
QD/nonlinear family, k-eigenvalue chapter, 14 symptom→cause→fix triples).

⚠ Notation: A&L normalize **Σ_n w_n = 2** (so ω_SI = (c/2)Σw_n/(1+μ_n²Λ²) → c at Λ=0,
Eq. 3.29, and their half-range γ_N = Σ_{μ_n>0}μ_nw_n → ½); Larsen 1982 and Alcouffe
1977 use **Σω = 1** (γ → ¼). Same objects, factor 2 — normalize before comparing rows.
Their λ is in MEAN-FREE-PATH units (modes e^{iΣ_tλx}).

### 4.1 The reference spectral results (the gate formulas)

- **SI**: ρ = c continuous (2.17, p. 26); discrete DD, ANY N and h (3.28–3.30, p. 51):
  ω_SI(Λ) = (c/2)Σ_n w_n/(1+μ_n²Λ²) with **Λ = (2/(Σ_th))·tan(Σ_thλ/2)** — the
  discrete map merely reparametrizes the mode axis; heterogeneous predictor
  min c ≤ ρ ≤ max c (3.31).
- **DSA continuous** (2.47–2.51, p. 29):
  ω(λ) = [3c/(λ²+3(1−c))]·[(λ²/3 + 1)(tan⁻¹λ)/λ − 1], σ_DSA ≤ **0.2247c**; at c = 1
  the sup sits at intermediate λ ≈ 2.5 (Fig. 1, p. 30); the c-scaling is sharp only
  for c ≈ 1. Iteration economics (2.52–2.54): ℓ_DSA ≤ 1.54m for 10⁻ᵐ; c = 0.999 ⟹
  ℓ_SI = 2301m.
- **Consistent discrete DSA (four-step), DD slab** (3.65, p. 57):
  ω = ω_SI − c(1−ω_SI)/(1 − c + ⅓Λ²) for ALL Σ_th — Σ_th-INDEPENDENT, and "for all
  quadrature orders N and all 0 < c ≤ 1.0 satisfy σ_DSA < 0.2247" (p. 57).
- **Inconsistent discrete DSA** (cell-edge conventional scheme; 3.41/3.43/3.44,
  p. 53): ω = ω_SI − c(1−ω_SI)/{(1−c)[1 + (Σ_th/2)²Λ²] + ⅓Λ²} — monotone-increasing
  in Σ_th; ω → ω_SI as Σ_th → ∞ ("convergent but ineffective — it performs like Source
  Iteration"); and for cell-CENTERED diffusion discretizations, outright divergence
  σ_DSA > 1 at large Σ_th [51 = Reed 1971] (p. 54).
- **Multi-D**: continuous results = 1-D exactly (plane-wave rotation, p. 81); discrete
  ordinates break rotational invariance ⟹ orientation-dependent ω. Consistent DSA
  (4.65–4.67, p. 83): σ = c/2 at S2, monotonically DECREASING to 0.2247c as N → ∞ —
  opposite to 1-D (0 at S2, increasing to 0.2247c). S₂SA multi-D: c/(2−c) (4.68).
- **Sphere**: continuous SI and DSA = slab EXACTLY (integral-equation trick,
  4.35–4.43, pp. 78–79; "this analytic result is new"); NO discrete curvilinear
  Fourier analysis exists (p. 79); cylinder: no theory at all (⅓ page, p. 80).

### 4.2 The consistency taxonomy (§III.B, §IV.D; Discussion pp. 138–139)

| Class | Definition | Stability | Examples |
|---|---|---|---|
| **Fully consistent** | low-order derived from the discrete SN equations by the discrete P1 projection (four-step) | unconditional; σ < 0.2247c (1-D); no known instability, ever (p. 69) | Alcouffe DD; Larsen four-step WD/LC/LD/LM; multi-D DD 9-point (2-D) / 27-point (3-D) vertex stencil — "yields a DSA spectral radius equal to that without spatial discretization" (p. 85) |
| **Partially consistent** | consistent derivation then simplified for solvability | works except on the hardest problems; degradation windows | TWODANT 5/7-point ("a winning strategy for problems that are not too optically thick and highly scattering," p. 87); M4S [169]; Wareing WDSA; asymptotic DSA (fails as cell aspect ratio → ∞, σ → c, p. 87); M4S diverges on skewed tets [300] (p. 88) |
| **Inconsistent** | independently-discretized diffusion operator | Σ_th-degradation (edge schemes) or divergence (cell-centered); "no safe dose" (Part II Table II) | Reed's scheme; Gelbard-Hageman 5-point (≡ edge scheme in slab) |

The canonical WHY (p. 58, verbatim): "In DSA, the importance of a consistently-
discretized low-order diffusion equation cannot be overstated. Generally, DSA with a
conventionally-discretized low-order diffusion equation becomes inefficient or even
divergent as Σ_th increases." Operator-level cause (p. 21): the preconditioner M must
have LARGE eigenvalues exactly where I − X_SI has small ones; overshoot ⟹ divergence —
"the operator M must be chosen very carefully." The engineering burden is "the tyranny
of consistent low-order discretization" (p. 20); the four-step = Lebedev's KP1(1)
applied to the fully-discrete problem, M = E D⁻¹K with K = discrete-P1 projection,
D = K(L−S)E the projected operator, E = prolongation (1.37–1.40, p. 58) — **the
moment-reduction statement of "consistent" that ORPHEUS's 3a computational derivation
targets.**

Escapes from the tyranny (each with a price; §II.C/E, §III.C, p. 138): S₂SA (same
spatial scheme for an S2 low-order — no consistency work, but ≡ DSA only in 1-D,
c/(2−c) in multi-D); TSA (low-order = damped transport, discretization-independent;
β-tradeoff, divergence risk at β → 1); nonlinear QD/CMFD family (no consistency
requirement for stability, but converged answer ≠ unaccelerated SN answer unless
consistent; positivity demands); second-order (even-parity/SAAF) forms (consistency
easy, but no sweep, banded solves, void trouble).

### 4.3 The reflective/periodic-boundary record (what the review actually contains)

The review does NOT develop an "inconsistent DSA at reflective boundaries" theory. What
it carries:

- The Fourier worst case is the INFINITE MEDIUM; in finite vacuum-bounded problems
  "the observed spectral radius is somewhat less than c … The inclusion of leakage
  shortens neutron lifetimes and hastens the convergence … this effect diminishes as
  the physical system becomes increasingly optically thick" (p. 51). Reflective/
  periodic problems eliminate leakage ⟹ they REALIZE the model-problem worst case —
  the operative reading for why acceleration trouble surfaces first on reflected
  problems. Lebedev's KP results were themselves derived on periodic model problems
  (p. 14).
- Yavuz & Larsen 1988, "Diffusion Synthetic Acceleration for S_N Problems with
  Reflecting Boundaries" (Trans. ANS 56, 305) is bibliography entry [135] but is never
  discussed in the body — the primary reflective-BC DSA source is OUTSIDE the review
  and outside the local library (flag for acquisition if the reflective-BC rows become
  contentious).
- Zika & Adams TSA with OPPOSING reflecting boundaries (p. 111, [284, 285]): the
  low-order problem then carries boundary angular-flux corrections in addition to
  scalar-flux corrections; the combined matrix is nonsymmetric; they recovered CG
  compatibility by identifying "the underlying symmetric problem (the infinite-medium
  problem simulated by the reflecting boundary conditions)."
- The exact discrete reflective rows for slab DSA are in Larsen Part I [(35): φ₁ = 0
  starting solve; (39a/b): f₁ = 0 correction solve] and were exercised throughout
  Part II (reflecting-left slab; Table II showed BC-inclusive inconsistency diverging).

### 4.4 DSA as a Krylov preconditioner (§VI, pp. 103–113; Discussion p. 140)

- Identification: SI ≡ Richardson on (I − L⁻¹S)ψ = L⁻¹q; DSA ≡ preconditioned
  Richardson, P = I + MS (1.27–1.34); scalar-flux operator T = I − X_SI, Tφ = η̂
  (6.4–6.5) with the matrix-free matvec Td = d − K₀L⁻¹(Σ_s/2)d (6.26) — exactly
  ORPHEUS's GMRES-on-sweep-preconditioned system (#200). Unifying identity:
  cond(T) κ ≈ 1/(1 − σ_SI) (6.20) — "finding a good preconditioner … is the same
  problem of finding a good low-order operator" (p. 106).
- **Krylov rescues inconsistent DSA** (Derstine-Gelbard Table 1, p. 110; XY, DD
  transport, cell-centered FD diffusion, c = 0.99): at Σ_tΔx = 0.5/1.0/2.0/4.0 —
  TWODANT (consistent) 6/6/7/12; inconsistent+CG 10/12/13/18; inconsistent alone
  12/15/div/div. "for some problems, inconsistent DSA diverges by itself but with CG
  it converges in less than 20 iterations!"
- **Krylov does not replace consistency** (Ashby et al., Table 2, p. 113; two-region
  slab, 50 mfp c = 0.99998 + 5×10⁷ mfp c = 0.8, 16-pt Gauss, 10⁻⁸, cap 50):

  | Transport | Diffusion precond. | Richardson | Chebyshev | GMRES |
  |---|---|---|---|---|
  | DD | consistent | 11 | 8 | 7 |
  | DD | FEM (inconsistent) | >50 | 22 | 13 |
  | LD | consistent | 11 | 8 | 7 |
  | LD | FEM (inconsistent) | 43 | 20 | 14 |

  Second Ashby problem (10⁸ mfp, c = 0.998): FEM preconditioner — ALL methods fail at
  50; consistent — all ≤ 2 (p. 113). Verdict: "nothing takes the place of a good
  preconditioner (such as a consistent diffusion discretization)" (p. 113); "a 'good'
  DSA scheme … converges so quickly, employing CG in addition generally reduces
  iteration counts by only one or two" (pp. 110–111).
- **The design guidance** (p. 140): use CG/GMRES *outside* DSA for problems with
  "severe spatial discontinuities, unstructured spatial grids, or even a slightly
  inconsistent acceleration scheme"; if the preconditioner suppresses all modes, the
  Krylov wrapper is "unnecessary — but probably not harmful." [The full quote is in
  the survey §2.7 — it licenses ORPHEUS's architecture: GMRES-with-sweep now (#200),
  DSA as the preconditioner upgrade (#2).]
- **Post-review thread** (the absent WWM 2004): the review already documents
  *consistent* multi-D DSA degrading to σ ≈ 0.88 on heterogeneous, high-aspect-ratio
  problems with "no known remedy" (pp. 69, 139 — degradation NOT instability;
  "(We are not aware that any problems have been seen in which a 'consistently-
  discretized' DSA scheme has actually become unstable.)" p. 69) and cites
  Warsa-Wareing-Morel as 2001 M&C proceedings [300]. The 2004 journal paper (Krylov
  restoring robustness there) is the natural continuation — absent locally, see §0.

### 4.5 Multigroup aspects (§I, §VII.A)

- What DSA accelerates: the WITHIN-GROUP self-scattering iteration — in ORPHEUS terms
  the ℓ = 0 (and ℓ = 1 under anisotropy) within-group moments. Group-to-group coupling
  is the outer (Gauss-Seidel over groups downscatter; power iteration for fission).
- Upscatter stagnation (§VII.A, pp. 114–117): zero-absorption thermal systems have a
  GS eigenvalue → 1 whose eigenvector is "a discrete approximation of a Maxwellian"
  (heavy water 41-group ρ_GS = 0.9998); fix = Adams-Morel Two-Grid [179, 180]:
  one-group diffusion correction weighted by the slow mode's spectral shape ξ_g
  (precomputed per material from the λ = 0 eigenproblem); ρ_TG = 0.489 (= the 2nd GS
  eigenvalue). Consistency carries over verbatim: "if one knows how to discretize DSA
  for a within-group scattering problem for a given transport discretization, one can
  apply the same discretization" to the TG equation (p. 117).
- k-eigenvalue DSA (§VIII.E): Alcouffe's nonlinear D̂/removal functionals = the
  DANTSYS production scheme; fixup interactions (M2 item 10); linear alternative
  (Gelbard et al. [81]) costs one adjoint low-order eigenvalue solve.

### 4.6 The diffusion-limit connection (why consistency, physically)

- The design principle (pp. 26–27): SI's slow modes are λ ≈ 0 — flat in space AND
  angle; the accelerator must be exact precisely there. A low-order operator whose
  discrete diffusion limit disagrees with the transport scheme's discrete diffusion
  limit leaves an O(1) residual on exactly those modes at finite Σ_th (Alcouffe's
  original R-does-not-vanish argument, §1.3 above — the review restates it as the
  (3.41) monotone degradation).
- Asymptotic DSA (Wareing-Larsen-Adams; §IV.D.4, pp. 86–87) makes the connection
  explicit by DERIVING the low-order operator from the transport scheme's asymptotic
  thick-diffusion-limit equations: effective in 1-D and on square 2-D cells; fails as
  aspect ratio → ∞ (σ → c) because "error modes that are discontinuous at cell edges;
  the diffusion-limit spatial discretization does not attenuate these modes" (p. 86).
  I.e., the diffusion-limit and the four-step routes coincide in 1-D but NOT in
  multi-D — consistency (moment reduction) is the stronger property.
- Related but distinct axis (NOT in this review's scope; see the local
  Larsen-Morel-Miller 1987 / Larsen-Morel 1989 papers + memory
  `space_angle_discretization_separability.md`): whether the transport scheme itself
  is ACCURATE in the thick diffusive regime. Do not conflate: DSA consistency is about
  the ITERATION converging; the diffusion-limit accuracy of the scheme is about the
  ANSWER being right.

### 4.7 Spatial-scheme specifics (DD first, LD arm gated)

- 1-D: the four-step is worked for WD/DD (review Eqs. 3.45–3.64 = Larsen Part I's DD
  instance), and Part I §IV–V supply LC/LD/LM. Differences from a conventional
  discretization (p. 56): a THREE-POINT removal term and the γ_N boundary factor.
  The four-step output depends on which algebraic form of the transport scheme you
  start from (which moments are taken) — p. 73 and the Fig. 8 p₀-vs-p₁ moment
  sensitivity (pp. 76–77): a warning that "the" consistent operator is
  presentation-dependent for non-DD schemes.
- Multi-D: DD eliminates to a 9-point (2-D) / 27-point (3-D) vertex diffusion stencil
  (p. 85); BLD/corner-balance schemes CANNOT be eliminated (12 unknowns/cell, p. 85) ⟹
  the modified-DSA family (M4S etc.) or Morel-Dendy-Wareing 1993 [190] multigrid-based
  BLD-consistent scheme ("converges any fixed-source one-group isotropic-scattering
  transport problem for the cost of a few transport sweeps," p. 88). First-order-form
  verdict (p. 69): the spatial discretization choice affects DSA convergence far more
  than the angular quadrature does.
- Unstructured grids: open as of 2002 (p. 89).

## 5. Morel 1982 — when φ0-only DSA stops being enough

J.E. Morel, NSE **82**(1), 34–46 (1982). Printed p.NN; PDF = printed − 32. Equations
verified: (16)–(22) on PDF pp. 4–5, (23)–(40) on pp. 5–8, (41)–(54) on pp. 8–10,
Appendix on pp. 13–14; Tables I–III on PDF pp. 11–12. Two contributions: (A) P1
acceleration (accelerate φ₀ AND J = φ₁) for anisotropic scattering; (B) "S_n diffusion
theory" (a curvilinear-corrected diffusion coefficient). Uses **Σw_m = 1**.

### 5.1 The operator framing (§II) — the cleanest general statement in the corpus

(A − B)ψ = Q; iterate Aψ^{l+1} − Bψ^l = Q; M = A⁻¹B; exact resummation
ψ = ψ^l + (A−B)⁻¹Br^l (Eq. 10, r^l = ψ^l − ψ^{l−1}); synthetic method = replace (A−B)⁻¹
by W⁻¹ (the low-order operator): ψ^{l+1} = ψ̃^l + W⁻¹Br̃^l (11b). Standard form (12):
ψ^{l+1} = [I − W⁻¹(A−B)]Mψ^l + (I + W⁻¹B)A⁻¹Q — "the spectral radius … goes to zero as
W gets closer to (A − B)." For transport: A = L (sweep), B = S, W = K = −∇·D∇ + σ_a
(15b). The scalar algorithm (16a–d): sweep; moments; φ^{l+1} = φ̃^l + K⁻¹[σ₀(φ̃^l −
φ^l)] (16c); and the **moment update** (16d): ψ^{l+1} = ψ̃^l·(φ^{l+1}/φ̃^l), realized
as: replace φ̃ by φ^{l+1}, SCALE all higher moments by φ^{l+1}/φ̃^l. Scaling "neither
increases nor decreases the convergence rate. It merely prevents the scattering source
from becoming negative" (p. 36). Eliminating: the corrected diffusion equation (20)

    −∇·D∇φ^{l+1} + σ_aφ^{l+1} = q − (∇·J̃^l + ∇·D∇φ̃^l)

— **independent confirmation of the resolved Alcouffe sign convention** (§1.5): defect
R = ∇·J̃ + ∇·D∇φ̃, SUBTRACTED. Alcouffe's three variants restated as (20)/(21)/(22).

**The sharpened consistency statement** (p. 37, load-bearing): "Consistency is a
necessary condition for convergence, but it is not a sufficient one. For instance, the
spatial differencing scheme used for the diffusion equation in early versions of the
synthetic method yielded a low-order equation that was CONSISTENT with the discrete
ordinates equations, but the resulting algorithm was nonetheless unstable. … [Larsen's]
work seems to indicate that unconditional stability can be expected if the low-order
equations can be DERIVED DIRECTLY from the S_n equations." I.e., fixed-point
compatibility (Reed had it) ≠ stability; the operative property is
derived-by-moment-reduction. Keep the two notions separate in the ORPHEUS docs.

### 5.2 P1 (anisotropic) acceleration (§III)

Low-order = the P1 operator on (φ, J) jointly (Eq. 23: off-diagonal ∇ and (1/3)∇;
diagonal σ_a and (σ_t − σ_1)). Algorithm (24a–c); update: even moments scaled by
φ^{l+1}/φ̃^l and odd by J^{l+1}/J̃^l (or all by the φ-ratio; identical convergence
rate, slightly better late-iteration balance for the even/odd split, p. 37).
P1-equivalent corrected diffusion (27) — Eq. (20) plus ONE extra correction component:

    −∇·D∇φ^{l+1} + σ_aφ^{l+1} = q − [∇·J̃^l + ∇·D∇φ̃^l + 3∇·σ₁D(J̃^l − J^l)] ,
    D = 1/[3(σ_t − σ_1)] ,

then the current update (28): J^{l+1} = J̃^l + D∇(φ̃^l − φ^{l+1}) + 3σ₁D(J̃^l − J^l).
"the extension from diffusion acceleration to P1 acceleration is both straightforward
and economical" (p. 38). Equivalent to Miller's generalized-rebalance derivation in
slab; **the curvilinear obstruction** (p. 38): generalized rebalance requires the SN
moments to satisfy the first two spherical-harmonic equations, but "in curvilinear
coordinates, the second equation is never satisfied because of the low-order angular
finite differencing scheme used to approximate the streaming term" (the α-recursion) —
the synthetic (moment-reduction) route has no such difficulty. [Direct relevance to
ORPHEUS's curvilinear SN: the four-step logic survives the α-recursion; rebalance-style
derivations do not.]

### 5.3 Discrete forms (slab + general 1-D)

- Slab discrete (29)–(34): diamond-differenced P1 low-order [(29a) balance + (29b)
  Fick with (29c/d) diamond averages] → cell-centered current (30) → two-cell sum (31)
  → the edge-flux diffusion system (32) [= Alcouffe's Eq. (21)-form × overall factor
  2, with ½σ_a(edge-sum)h removal]. Corrected versions (33)/(34a) add the anisotropic
  term −6[σ₁D_{i+1}(J̃_{i+1} − J_{i+1}) − 6σ₁D_i(J̃_i − J_i)]-type differences; current
  update (34b).
- **General 1-D geometry** (35)–(38): with areas A_{i±1/2} and volumes V_i, summing two
  adjacent balance rows does NOT eliminate the shared edge current — remainder
  (A_{i+3/2} − A_{i−1/2})J_{i+1/2} (Eq. 36). Following Alcouffe, approximate the
  balance with cell-centered areas A_i = ½(A_{i+1/2} + A_{i−1/2}) (37a/37b) —
  "approximate ONLY for spherical and cylindrical geometry" (exact in slab) — giving
  the curvilinear corrected diffusion (38).
- Fixup consistency (39)/(40): either the LINEAR extra correction (39) (≡ 0 under pure
  diamond) or **Alcouffe's nonlinear removal** (40a/40b):
  σ̃_a^l = (σ_aφ̃^l_{i+1}V_{i+1} + σ_aφ̃^l_iV_i)/φ̃^l_{i+1/2}. Morel USES the nonlinear
  form to stay equivalent to Alcouffe's source-correction scheme under isotropy.
- **Stability caveat** (pp. 39–40, verbatim-critical): "neither general geometry,
  negative flux fixup consistent scheme is unconditionally stable. Alcouffe's proof of
  stability rigorously applies only to Eqs. (34a) and (34b) with σ₁ = 0. It does not
  apply to Eqs. (37a) and (37b) even for the diamond-difference slab geometry case
  because of the nonlinear absorption term." Empirically stable absent significant
  fixup activity. ⟹ The UNCONDITIONAL-stability guarantees in this corpus attach to
  the LINEAR slab schemes only; curvilinear and nonlinear-removal variants are
  empirically-supported, not proven.

### 5.4 Boundary rows (Appendix, pp. 45–46) — the only explicit ALBEDO rows in the corpus

Boundary-cell equations for the discrete corrected diffusion (40a), left cell (right
analogous), one row covering vacuum/reflective/albedo via α (α = 0 vacuum, α = 1
reflective, α = albedo value):

    2[−A₁(D₁/h₁)(φ^{l+1}_{3/2} − φ^{l+1}_{1/2}) + 0.5·A_{1/2}·((1−α)/(1+α))·φ^{l+1}_{1/2}]
    + σ̃_a* φ^{l+1}_{1/2}
    = q₁V₁ − (A_{3/2}J̃^l_{3/2} − A_{1/2}J̃^l_{1/2})
      + 2[−A₁(D₁/h₁)(φ̃^l_{3/2} − φ̃^l_{1/2}) + 0.5·A_{1/2}·((1−α)/(1+α))·φ̃^l_{1/2}]
      − 6A₁σ₁D₁(J̃^l₁ − J^l₁) ,       σ̃_a* = σ_aφ̃^l₁V₁/φ̃^l_{1/2} .

The (1−α)/(1+α) factor is the Marshak-type albedo current coefficient (vacuum ½·φ
term; reflective 0). A periodic row is also given (couples cells 1 and K). Note
(p. 46): "the vacuum condition is used for the diffusion equation at a vacuum boundary
EVEN IF boundary sources are present in the discrete ordinates calculation" (boundary
sources enter the transport sweep only; optionally mock them as a first-cell
distributed source for the uncorrected initial-guess solve). ⚠ The correction terms on
the RHS mirror the LHS operator applied to φ̃ — same defect-cancellation structure as
the interior rows; and the OCR of the last-cell row shows index slips (A_{K−1/2} vs
A_K, σ_K vs σ₁) — re-verify against PDF p. 14 before transcribing THAT row into code.

### 5.5 S_n diffusion theory (§IV) — curvilinear D correction

Derived from the curvilinear discrete SN equations (41) with the standard α-recursion
(42) [α_{m+1/2} = α_{m−1/2} − μ_mw_m, α_{1/2} = α_{N+1/2} = 0], STEP differencing in
angle (44), and the linear-in-μ ansatz ψ_m = φ + 3Jμ_m (45). The P1 moment yields
Fick's law (50a) with the modified coefficient (50b):

    D_i = A_i h_i / [3(A_{i+1/2} − A_{i−1/2})β + 3(σ_t − σ_1)V_i] ,
    β = 3Σ_m α_{m−1/2}(μ_m² − μ_mμ_{m−1})   (49b) ,

h → 0 limits (51a–c): D = 1/[3σ_tr] (slab); 1/[3(σ_tr + β/r)] (cylinder);
1/[3(σ_tr + 2β/r)] (sphere) — D → 0 as r → 0. Physical content: near a point source
the SN transport solution has an r⁻² component; SN-diffusion reproduces it (53b) while
standard diffusion does not (53a). With angular DIAMOND differencing β ≈ 0 (theory ≈
standard); β → 0 as N → ∞ anyway. Practical value (Table I): mild — point-source
sphere S4: 8 → 5 iterations; distributed source: 5 → 4; N = 16: no difference. Morel's
caution (p. 45): this claims nothing about SN-diffusion being more ACCURATE generally
— it is a convergence-acceleration device.

### 5.6 Measured P0-vs-P1 envelopes (§V; Tables II–III, PDF pp. 11–12)

Setup: NONABSORPTIVE medium (c = 1 — worst case), vacuum boundaries, 30 uniform cells,
Gauss quadrature, 10⁻⁴ pointwise, S4 unless stated.

**Table II** (P1 scattering, cell = 1 TRANSPORT mfp, so σ₀h = 1/(1−μ̄) grows with
anisotropy; "isotropic" = P0-DSA, "anisotropic" = P1-DSA):

| μ̄ = σ₁/σ₀ | σ₀h | slab P0 / P1 | sphere P0 / P1 |
|---|---|---|---|
| 0.0 | 1.0 | 4 / — | 4 / — |
| 0.25 | 1.33 | 8 / 4 | 10 / 4 |
| 0.50 | 2.0 | 14 / 4 | 15 / 4 |
| 0.75 | 4.0 | 31 / 4 | 31 / 4 |

P1-DSA count is FLAT at 4; P0-DSA degrades steadily with μ̄. "the effectiveness of P1
acceleration with P1 scattering is equivalent to that of P0 acceleration with P0
scattering" (p. 43).

**Table III** (Fokker-Planck expansion of degree N−1 with S_N — μ̄ → 1 with N; cell =
⅓ transport mfp; sphere runs use S_n diffusion theory; slab: isotropic boundary
source; sphere: point source):

| σ₀h | N | μ̄ | slab P0 / P1 | sphere P0 / P1 |
|---|---|---|---|---|
| 2.0 | 4 | 0.8333333 | 42 / 9 | 42 / 10 |
| 9.33 | 8 | 0.9642857 | 142 / 26 | 144 / 40 |
| 40.0 | 16 | 0.9916666 | 374 / 72 | 482 / 127 |

P1-DSA beats P0-DSA by 4–6×, but ALSO degrades as the ℓ > 1 moments strengthen
(72–127 at S16) — the ladder does not end at P1. Operational incidents (p. 44): the
sphere S8/S16 runs needed spatial STEP differencing; the ONEDANT fixup interacted with
the acceleration in an OSCILLATORY manner (fixup toggling on/off per iteration) — with
diamond, S8 took 80/166 iterations, S16 failed at 1000; oscillatory fixup
destabilization is "greatly increased in problems in curvilinear coordinates having
cell thickness on the order of ten total mean-free-paths or more" and was also seen
with ONETRAN's LD + rebalance (not specific to DSA).

### 5.7 The "is φ0-only DSA sufficient?" criterion for ORPHEUS

- Morel's paper gives the trend (Tables II/III) and the design rule; the QUANTITATIVE
  mode-by-mode criterion is A&L 2002 Eq. (7.14) (p. 119): at λ = 0 the per-moment SI
  eigenvalues are ρ_n = Σ_sn/Σ_t; DSA removes the n = 0 (P0-DSA) or n ≤ 1 (P1-DSA)
  modes, leaving ρ ≈ max_{n>accelerated} Σ_sn/Σ_t. So: **P0-only DSA is sufficient when
  Σ_s1/Σ_t is comfortably below the target spectral radius (~0.23); switch to P1
  acceleration when Σ_s1/Σ_t approaches Σ_s0/Σ_t; expect degradation of even P1-DSA as
  Σ_s2 → Σ_s0** (McCoy-Larsen §IV restatement; angular-multigrid territory beyond
  that, A&L §VII.B).
- ORPHEUS carries Legendre-anisotropic scattering ⟹ the DD-DSA arm should implement
  P0-DSA first (matches the isotropic-scattering unconditional-stability theory), with
  the P1 extension designed-for (Larsen Part I already carries φ₁ acceleration through
  every scheme; Morel (27)/(28) is the compact continuous statement) and gated on
  measured Σ_s1/Σ_t in the target problems.

## 6. Synthesis — mapping onto the ORPHEUS Phase-3 build

Maps the extracted math onto the build plan WITHOUT designing the implementation.

### 6.1 "Consistent with DD" as a discrete-matrix statement

Supported by the sources at three levels, strongest last:

1. **Fixed-point compatibility** (necessary only): the low-order equation's converged
   solution equals the discrete SN solution — Alcouffe Eqs. (19)/(20); Morel p. 37.
   Reed's unstable scheme HAD this. Not the operative property.
2. **Diffusion-limit compatibility** (Alcouffe's diagnostic): the correction source
   R (the diffusion defect of the transport iterate) must vanish IDENTICALLY on the
   discrete diffusion limit of the SN scheme at ANY h — realized by R = second
   difference of φ̃₂ (Alcouffe Eq. 27). Cell-centered differencing fails this because
   the DD P1-moment equation (16) staggers current (centers) against scalar flux
   (edges).
3. **Moment reduction / Schur complement** (the operative definition; the target of
   ORPHEUS's 3a): the low-order operator is DERIVED from the assembled discrete SN
   operator by the discrete P1 projection — A&L (1.37)–(1.40): M = E D⁻¹ K with
   K = the discrete-P1 analysis (take L₀ AND L₁ moments of BOTH the balance and the
   closure — Larsen (16a–d)), E = the linear-in-angle prolongation
   (ψ = φ₀ + 3μφ₁, Larsen (15a)/(33)), D = the projected operator "K applied to
   L − S". Concretely for slab: form the four moment equations; assign iterate
   indices per Larsen step 3 (l+1 on all φ₀/φ₁ slots; lag φ₂ and the closure
   residuals L₀γψ, L₀βψ at l+1/2); subtract the un-promoted system; **Schur-eliminate
   the current corrections f₁ and the cell-average corrections** onto the edge
   scalar-flux corrections f₀ ⟹ the tridiagonal Larsen (27) with coefficients (23a–f)
   [DD: ρ_i = 0 ⟹ D_i = 1/(3σ_T), σ̂_R = σ_R, three-point ¼-weighted removal]. The
   boundary rows come from the SAME reduction applied to the half-range (Marshak)
   moments: (38a/b) with the discrete γ_N = Σ_{μ>0}μ_mω_m, (39) reflective f₁ = 0,
   (40) periodic, closed via the one-sided relations (25)/(26).

   **The correspondence check for 3a**: the computationally-derived low-order matrix
   must equal Larsen (27) + (28) + BC rows entry-for-entry (DD instance), and its
   convergence identity must be the literal "0 = 0" collapse (Larsen p. 54). Watch
   items: the ⅓ in the P1 moment (OCR trap, §1.2); the sign of the defect source
   (Alcouffe print slips, §1.5); the ¼ three-point removal weights; γ_N ≠ ¼
   (Σω = 1) at small N; the residual source g₀ = σ̂_S h (φ₀^{l+1/2} − φ₀^l).

4. Presentation-dependence caveat (A&L p. 73): for non-DD schemes the four-step
   output depends on WHICH algebraic form of the transport scheme is reduced (p₀ vs
   p₁ moments, Fig. 8). For DD the result is canonical. Record the chosen form when
   the LD arm lands.

### 6.2 The reference recipe for slab DD (what 3a must reproduce)

Primary: Larsen Part I §II–III [(10)→(16)→(18)→(19)–(27)→(28) + (30)–(40)] at
ρ_i = 0. Secondary (algebraically equivalent, nonlinear-removal form): Alcouffe
(22)/(23)/(23a) with the §1.5 sign resolution; under the diamond assumption the
(23a) removal linearizes to Larsen's ¼-weights. Cross-check: A&L (3.45)–(3.64)
(their Σw = 2). Morel (32)/(34) is the same operator × 2 with the P1 extension.

### 6.3 The ρ gates (formula, value, validity)

| Gate | Prediction | Source | Valid when |
|---|---|---|---|
| SI | ρ = c; discrete ω_SI(Λ) = (c/2)Σw_n/(1+μ_n²Λ²), Λ = (2/Σ_th)tan(Σ_thλ/2), Σw = 2 | A&L (2.17), (3.28)–(3.30) | inf. hom. medium, isotropic, any N/h; finite media: ρ ≲ c (leakage helps); heterogeneous: min c ≤ ρ ≤ max c (3.31) |
| SI+DSA continuous | ρ ≤ 0.2247c; ω(λ) closed form (2.50) | A&L (2.50)/(2.51) | continuous limit (fine mesh), isotropic, P1 low-order |
| SI+DSA consistent DD | ω(Λ) = ω_SI − c(1−ω_SI)/(1−c+⅓Λ²); σ < 0.2247, Σ_th-INDEPENDENT | A&L (3.65) | inf. hom., isotropic, uniform h, four-step DD operator — **the primary quantitative gate** |
| WD/LC/LD/LM model bound | spr ≤ cK_N; Table I (WD: 0.185→0.248; LC ≤ 0.315; LD ≤ 0.300; LM ≤ 0.282; K_N = 1−3ρ² at the WD sup; K_∞ = ¼) | Larsen Tab. I, p. 51 | model problem; sup over the α-family — pure DD at finite h sits BELOW (use 3.65) |
| S2 exactness | K_2 = 0: converge in ONE iteration (to roundoff) | Larsen p. 56, McCoy-Larsen (1b) | requires exactly consistent BC rows — the sharpest BC unit test |
| Inconsistent control | ω = ω_SI − c(1−ω_SI)/{(1−c)[1+(Σ_th/2)²Λ²]+⅓Λ²} (edge scheme, degradation); cell-centered: divergence at large Σ_th | A&L (3.41)/(3.43)/(3.44); Reed | negative-control: the battery must SEE this |
| Iteration counts | ≤ ~8 iters to 10⁻⁴ pointwise at ANY Σ_th (c = 0.98, S4); "two orders of magnitude per three iterations" | McCoy-Larsen Tab. I/III; A&L p. 58 | consistent scheme, no fixups |
| ρ measurement | σ̂ ≈ ‖φ^{l+1}−φ^l‖/‖φ^l−φ^{l−1}‖; stop-test ε(1−σ̂) | A&L (1.15)–(1.20) | dominant eigenvalue real+simple (usual) |

### 6.4 Stability discriminators for the test battery

1. **Cell-thickness sweep**: σ_th ∈ {1, 2, 4, 6, 10, 20} (McCoy-Larsen grid; Alcouffe
   2-D went to 15 mfp; A&L Table 2's second problem to 10⁸ mfp total). Consistent DSA:
   flat counts; inconsistent: monotone degradation/divergence.
2. **c → 1**: c ∈ {0.8, 0.98, 1.0} (pure scatterer exercised by Alcouffe Table III and
   Morel throughout); plus the c = 0.99998/50-mfp Ashby-style two-region slab (A&L
   Table 2 spec) as the hard case.
3. **Reflective BC**: reflecting-left slab (McCoy-Larsen standard problem);
   reflective realizes the no-leakage worst case (A&L p. 51 argument). Include a
   both-sides-reflective variant with c < 1 (note: c = 1 + all-reflective is
   physically ill-posed for fixed-source — the correction operator's (1−c)σ removal
   is what keeps it invertible; watch conditioning as c → 1).
4. **Partial-consistency negative control**: WD transport (a = 0.5) + DD-derived
   low-order INCLUDING its BC rows — expect Table-II-shaped divergence at σ_th ≳ 10
   (McCoy-Larsen). Proves the battery detects inconsistency (both operator AND BC).
5. **Anisotropy ladder**: μ̄ = σ_S1/σ_S0 ∈ {0, 0.25, 0.5, 0.75} at fixed transport-mfp
   cell width (Morel Table II): P0-DSA counts should track ~Σ_s1/Σ_t degradation;
   flat if/when P1-DSA lands.
6. **Heterogeneous regression**: the McCoy-Larsen shielding problem (§3.1 — full spec
   transcribed) incl. the coarsest (1,1,2,1) mesh; and Alcouffe's Table-III 2-D square
   (pure scatterer, reflective×2) when multi-D arrives.
7. **Fixup interactions**: N/A if ORPHEUS SN carries no negative-flux fixup (the
   correctness-first posture); if any fixup/limiter exists, McCoy-Larsen §III is the
   hazard map (fixup on ≈-negative meshes = "anathema"; oscillatory toggling in
   curvilinear ≥ ~10 mfp cells, Morel p. 44).

### 6.5 Notation/normalization mapping (sources → ORPHEUS)

- Slab μ_m = ORPHEUS `mu_x`; edge/center staggering i±1/2 vs i as in the SN mesh.
- Quadrature normalization: Alcouffe/Larsen/McCoy/Morel **Σw = 1** (φ = average;
  ⟨μ²⟩ = ⅓; half-range γ_N → ¼; Marshak ½φ₁ coefficient); A&L 2002 **Σw = 2**
  (γ_N → ½). Map to ORPHEUS's quadrature convention ONCE, at the start of 3a, and
  assert Σw and Σwμ² numerically in the derivation harness.
- Iterate indices: l+1/2 = post-sweep (Alcouffe's tilde), l+1 = post-low-order-solve
  — matches A&L and the ORPHEUS DSA plan vocabulary.
- σ_R = removal = σ_t − σ_s0,g→g (within-group); c = σ_s/σ_t; D = 1/(3σ_tr).
- The f-form (correction) vs full-flux form: algebraically equivalent (Larsen §II);
  prefer the f-form (storage, "0 = 0" convergence identity, fixup-insensitivity,
  homogeneous BCs). Alcouffe's production form solves for the full flux with the
  defect-corrected source and a NONLINEAR removal — historical, not the recommended
  target.
- Operator vocabulary: the sources' L includes collision (their L = ORPHEUS L + C);
  their sweep-preconditioned A = I − L⁻¹S and T = I − X_SI are NOT the ORPHEUS honest
  matvec A = L + C − S − B (see the A&L survey's notation block).

### 6.6 Open questions for the plan-of-record

1. **Reflective-BC low-order rows**: Part I gives the exact slab rows (f₁ = 0 closed
   via (25)/(26)); no source ANALYZES reflective-BC DSA stability (Yavuz-Larsen 1988
   [135] is the missing primary — Trans. ANS summary, absent locally). Empirics
   (McCoy-Larsen reflecting-left; Alcouffe 2-D reflective×2) show no anomaly for
   consistent schemes.
2. **Near-singularity at c → 1 with all-reflective/periodic BCs**: the correction
   operator's removal ∝ (1−c) ⟹ conditioning grows; none of the sources treat it
   (their reflected problems keep a vacuum side or c < 1). Decide: guard, regularize,
   or exclude (fixed-source c = 1 no-leakage is ill-posed anyway).
3. **Eigenvalue coupling** *(rev. 2, Adams-Martin/WWM)*: Part I/II are fixed-source;
   Alcouffe's k-eigenvalue variants are nonlinear (D̂/removal correction, DANTSYS
   lineage) with fixup interactions; the linear alternative (Gelbard et al. [81],
   A&L §VIII.E) costs an adjoint low-order eigenproblem — and Adams-Martin's Ref. 19
   (Adams & Larsen 1989 Santa Fe proc., "Linear Diffusion-Synthetic Acceleration of
   k-Eigenvalue Problems") is now a NAMED primary for that linear-DSA-k route (§7.6).
   **The third route is now extracted firsthand (§8)**: keep the k-loop
   linear-algebraic — Arnoldi/IRAM on A = DH⁻¹F (moment-vector eigenproblem,
   parameter-free, outer counts flat in the dominance ratio: 14–26 outers across
   δ = 0.970–0.999 vs 468–8824 for power iteration, Warsa et al. Table IV) — and
   leave DSA entirely at the within-group level, in either posture (SI+DSA inners
   for a power-iteration loop; DSA-preconditioned FGMRES/BiCGStab inners + the BFG
   relaxed-tolerance cascade for an Arnoldi loop, which forfeits warm starts).
   ORPHEUS #2 scope = within-group/fixed-source acceleration first; the k-loop
   design is a separate decision — but it is now a DESIGNED-FOR decision with a
   production-scale precedent (28M DOF; upscatter nesting is the recorded weak spot).
4. **Multigroup/outers**: within-group ℓ = 0 DSA is the Phase-3 core; outer/upscatter
   acceleration (Alcouffe Eq. 5 implicit-multigroup diffusion; Adams-Morel Two-Grid)
   is future work — consistency discretization carries over verbatim (A&L p. 117).
5. **Curvilinear**: NO discrete curvilinear Fourier theory exists (A&L p. 79); the
   practice is Morel's cell-centered-area approximation (37a) with NO unconditional
   proof (Morel pp. 39–40) + empirical stability. The ORPHEUS curvilinear arm will be
   empirically gated only; the α-recursion blocks rebalance-style derivations but not
   the moment-reduction route (Morel p. 38).
6. **Multi-D consistency cost** (when 2-D arrives): DD eliminates to 9-point vertex
   (Alcouffe §II.C corner moments); non-DD needs the modified-DSA family or
   Morel-Dendy-Wareing multigrid; production codes accept slight inconsistency
   (TWODANT 5-point) + optional Krylov wrap (A&L p. 140).
7. **Absent primaries** *(rev. 2, Adams-Martin/WWM)*: ~~Adams & Martin 1992~~ now
   LOCAL + extracted (§7); the k-eigenvalue Warsa et al. NSE 147(1) also LOCAL +
   extracted (§8). Still absent: **Warsa-Wareing-Morel NSE 147(2), 218 (degraded
   consistent DSA + Krylov fix — the §0 disambiguation applies)** and **Yavuz-Larsen
   1988 (reflective-BC DSA)**. User to acquire if/when those arms open.
8. **LD-arm reference** (gated) *(rev. 2, Adams-Martin/WWM)*: TWO published routes —
   (a) fully consistent: Larsen Part I §V — LD α_m = ε_m/(3+|ε_m|), tridiagonal (78)
   with (81) coefficients, K_N ≤ 0.300 (coarse-mesh LD limit 0.25 per A&M printed
   156); plus the A&L p. 73 presentation-dependence warning; (b) M4S (§7): lag the
   first-moment interface-jump term C in step 2 ⟹ exact discrete Fick collapse
   (32), always-eliminable weak form (34)/(35), tridiagonal via App. C — simpler,
   dimension- and curvilinear-robust (spherical App. A), ρ → 0 in thick cells, at
   the cost of a 0.50–0.67c intermediate-thickness peak and the A&L-recorded
   skewed-tet divergence [300]. A&M verified the θ = 1 LD variant of Larsen-Morel
   1989 JCP 83 (the ORPHEUS slab-LD reference) behaves identically (§7.4). Choose
   when the arm opens; the §7.4 table is the decision surface.

### 6.7 Reading order for the implementing session

A&L §§II.B + III.B (tutorial + consistency mechanics, with the survey file) →
Larsen Part I §§II–III (the recipe + BCs — the 3a reference) → McCoy-Larsen (the gate
tables) → Alcouffe (DD specifics, 2-D corner scheme, eigenvalue variants; mind §1.5
sign errata) → Morel (P1 extension + curvilinear + albedo rows) — matching the survey
Task-5 recommendation.

---

## 7. Adams & Martin 1992 — the M4S paper (deliberately inconsistent DFE diffusion)

M.L. Adams & W.R. Martin, "Diffusion Synthetic Acceleration of Discontinuous Finite
Element Transport Iterations," NSE **111**(2), 145–167 (1992), DOI
10.13182/NSE92-A23930. LOCAL. Page map: printed ≈ PDF + 143 (printed 145 = PDF 2).

> **Provenance note (rev. 2)**: the Mistral OCR key returned persistent 401s on
> 2026-07-26 (worked at 00:50 the same day — key expired/rotated mid-day), so this
> paper's search surface is a locally-generated pypdf text-layer file
> `scratch/literature_ocr/<stem>.textlayer.md` (same `## p. N` = PDF-page convention),
> NOT a Mistral sidecar. All load-bearing equations below were verified against the
> rendered PDF pages 8–13 and 15–17 (`Read pages=`). Regenerate the real sidecar once
> the key is fixed (`tools/ocr_literature.py --glob 'Adams-Martin*'` — cache-idempotent).

### 7.1 What the paper is (abstract + §I, printed 145–146, PDF 2–3)

A discontinuous-finite-element (DFE) discretization of the diffusion equation for
accelerating DFE transport (linear discontinuous LD and bilinear discontinuous are the
worked instances), obtainable "in any coordinate system on any grid." The paper's own
headline claims (abstract, verified PDF 2): the diffusion discretization is **"not
strictly consistent with the transport discretization in the usual sense"**, yet the
scheme is **unconditionally stable and rapidly convergent**, and **ρ → 0 as cells
become optically thick** (instant convergence — unique among DSA schemes of the time;
shared only by boundary-projection acceleration, their Refs. 12/16/17, printed 160).
Analysis: 1-D + 2-D Cartesian. Numerics: 1-D Cartesian AND spherical. Motivation
(§II, printed 148, PDF 5): standard linear DSA (= Larsen's fully-consistent four-step,
their Ref. 2 = memo §2) is "algebraically complicated and potentially difficult to
solve … it may not be possible to eliminate the current vectors … in general
multidimensional problems. This complexity has precluded any widespread use."

### 7.2 The modified four-step — step 2 is the modified step (§IV.B + Table I, printed 152–153, PDF 9–10)

Table I (printed 152, PDF 9, verified verbatim) compares "Standard Linear DSA
Procedure" vs "Our Modified Four-Step Procedure":

| Step | Standard linear DSA | M4S |
|---|---|---|
| 1a | zeroth angular moment of the discretized transport eq. | Same |
| 1b | first angular moment | Same |
| 2 | promote indices to l+1 **except on second-and-higher moment terms** | "Same, except do not change indices on certain zeroth- and first-moment terms." |
| 3 | subtract acceleration eqs from unaccelerated eqs | Same |
| 4 | eliminate first moments → discretized diffusion eq. "May not be possible given high-order discretization schemes in two or three dimensions." | "Same; always possible" |

**The one modification**: after steps 1a/1b produce the discrete moment equations
(20) [zeroth] and (27) [first] — *identical to standard linear DSA up to this point*
(stated twice, printed 153) — step 2 ALSO leaves at l+½ the **surface term C(r_k)**
of the first-moment equation (27). Verbatim (printed 153, PDF 10): "Our procedure is
the same except that we leave the index at (l+½) on the term C(r_k) that appears in
Eq. (27). Because this term contains contributions from scalar fluxes and currents,
our scheme is not consistent."

The lagged object (unnumbered display after Eq. (27), printed 153, PDF 10, verified
on the render):

C^{(l+½)}(r_k) = (n_k/6)[φ₀(r_k⁺) − φ₀(r_k⁻)] + (3/16)[**φ₁**(r_k⁺) − **φ₁**(r_k⁻)]
+ (3n_k/16) n_k·[**φ₁**(r_k⁺) − **φ₁**(r_k⁻)] + h.o.t.^{(l+½)}

— i.e. the DFE **interface-jump** (discontinuity) contributions of scalar flux and
current in the first-moment equation. The first-moment volume operator itself is
promoted normally: (30) carries ⅔∇·φ₂^{(l+½)} (lagged, as standard) + ⅓∇φ₀^{(l+1)} +
σ_t φ₁^{(l+1)} — the ⅓ P1 coefficient renders correctly here (the §1.2 OCR trap does
not bite in this paper's print).

**Why this modification buys solvability**: with C lagged, EVERY surface term of the
first-moment equation subtracts out in step 3, and the volume terms collapse to an
exact discrete Fick law — Eq. (32), printed 154: **f₁(r) = −D∇f₀, D = 1/(3σ_t)**
[stated caveat: "we have assumed that gradients of our DFE basis functions reside in
the function space spanned by the basis functions. This is true of most DFE
methods."]. Hence step 4 (current elimination) is ALWAYS possible, in any dimension,
on any grid — the exact point where Larsen's fully-consistent reduction can jam
(Table I step-4 row; A&L p. 73 presentation-dependence is the same pain point).

**What "consistent" means here** (printed 152, verified): "The basic idea is to
derive a diffusion discretization that is equivalent to the transport discretization
when the transport angular flux is linearly anisotropic. This is what is meant by a
'consistent' diffusion discretization. Standard linear DSA achieves this rigorously,
while our procedure only approximates it."

**Vocabulary caution**: the paper NEVER says "partially consistent" — its own words
are "not strictly consistent" and "inconsistent, but rapidly convergent" (printed
153/160). The "partially consistent" bucket name the memo's §4.2 table uses is
A&L 2002's later taxonomy (their §IV.D, citing this paper as [169]). Credit the
taxonomy to A&L, the existence proof to this paper.

### 7.3 The resulting low-order operator (§IV.A ≡ §IV.B; Eqs. (15)–(19) and (34)–(36), printed 151/154–155, PDF 8/11–12)

**Two-path equivalence** (the paper's structural centerpiece, Fig. 1): discretizing
the CONTINUOUS diffusion equation by the same DFE machinery as the transport scheme
(same weights v_ki, same basis b_kj, upstream-evaluated partial currents) — Eq. (18)
— produces IDENTICALLY the operator that the modified four-step extracts from the
DFE transport discretization — Eq. (34). Verified statement under (34), printed 154:
"which is identical to the DFE diffusion equation that was derived by a completely
different procedure in Sec. IV.A."

The general weak form (Eqs. (18) ≡ (34), any grid, Cartesian; verified PDF 8/11):

∮_{∂Z_k} v_ki [(αf₀ − ½ n·D∇f₀)_{r_k⁻} − (αf₀ + ½ n·D∇f₀)_{r_k⁺}]
+ ∫_{Z_k} [D∇f₀·∇v_ki + σ_a v_ki f₀]
= ∫_{Z_k} v_ki σ_s [φ₀^{(l+½)} − φ₀^{(l)}],  1 ≤ i ≤ J_k, 1 ≤ k ≤ K

with the net→partial-current replacement (16): outgoing partial current evaluated
just INSIDE the zone (r_k⁻), incoming just OUTSIDE (r_k⁺) — "upstream" for partial
currents where it is undefined for the net current. Quadrature constants (verified
PDF 9): **(25a) Σ_{n·Ω_m>0} w_m|n·Ω_m|² = 1/6 and (25b) α ≡ Σ_{n·Ω_m>0} w_m|n·Ω_m|
≈ ¼** ⟹ **Σw = 1 normalization** (the Alcouffe/Larsen/McCoy/Morel camp of §6.5, NOT
A&L's Σw = 2). α is the discrete half-range current weight — the same object as
Larsen's discrete γ_N (memo §2.4), quadrature-exact rather than the continuous ¼.

**Boundary conditions** (19)/(26), verified PDF 8/10: incoming partial current on the
problem boundary = β·(incident-flux term) + (1−β)·(outgoing partial current); β = 0
reflective, β = 1 vacuum/specified-incident. In the correction form the incident term
cancels: (33b) ΔJ_in = γ(r_k)·ΔJ_out with γ = 1 − β (follows from differencing (26)).
⚠ Direction flag: A&M's β is an ABSORPTION weight (β = 0 ⟹ reflective), the reverse
sense of Morel's albedo α (memo §5.4, factor (1−α)/(1+α), α = 1 ⟹ reflective).

**Slab LD instance** (§IV.C, Eqs. (35a/b/c) + update (36), printed 154–155, PDF
11–12, verified): two unknowns per cell — cell-average f₀k and slope moment f₀k^x:

- (35a) balance row: −[D df₀/dx]_{k+½} + [D df₀/dx]_{k−½} + σ_ak h_k f₀k
  = σ_sk h_k [φ₀k^{(l+½)} − φ₀k^{(l)}]
- (35b) slope row: −3[D df₀/dx]_{k+½} − 3[D df₀/dx]_{k−½} + 12 D_k f₀k^x/h_k
  + σ_ak h_k f₀k^x = σ_sk h_k [φ₀k^{x(l+½)} − φ₀k^{x(l)}]
- (35c) edge net current from upstream partial currents, interior edge k+½:
  [α(f₀k + f₀k^x) − D_k f₀k^x/h_k] − [α(f₀,k+1 − f₀,k+1^x) + D_{k+1} f₀,k+1^x/h_{k+1}]
  (left boundary k = 0: first bracket → 0, incident BC; right k = K: whole
  expression → 0, reflective). The interior D_k f₀k^x/h_k terms are ½·D·(df₀/dx)
  with df₀/dx = 2f₀k^x/h_k for the LD basis — self-consistent with (16).
- (36): the update adds the correction to **BOTH** spatial moments:
  φ₀k^{(l+1)} = φ₀k^{(l+½)} + f₀k^{(l+1)} AND φ₀k^{x(l+1)} = φ₀k^{x(l+½)} + f₀k^{x(l+1)}.

**Solvability**: Appendix C manipulates (35) into three-point (tridiagonal) form in
EITHER cell-average unknowns ((C.3)–(C.5)) or cell-edge unknowns ((C.6)–(C.12)) — so
the 1-D M4S low-order solve costs the same as the memo §2's Larsen tridiagonal. In
2-D this reduction does NOT exist ("we cannot produce a difference equation involving
only one unknown per cell or vertex," printed 160) — the efficient 2-D solve was an
OPEN question, pointed at Morel's multigrid (Ref. 18, private comm. 1991 — the thread
that became Morel-Dendy-Wareing 1993 [190], already in memo §4.7/§6.6.6).

**Khalil distinction** (printed 151, verified): Khalil (NSE 98:226, 1988) also
discretizes the continuous diffusion equation mimicking the transport scheme, but
mimics its *balance + auxiliary cell relations*; A&M mimic its *weight + basis
functions*. Different diffusion operator, different DSA (compared in §V).

### 7.4 Reconciliation with the memo's LD reference (§2.6 Larsen Part I §IV–V)

**Verdict: a DIFFERENT operator, not a rederivation and not a simplification of the
same algebra.** Three low-order operators now sit on record for slab LD transport:

| Operator | Derivation | ρ fine-mesh | ρ intermediate | ρ thick limit |
|---|---|---|---|---|
| Larsen Part I (78)+(81) ("standard linear DSA", fully consistent) | four-step, unmodified | 0.2247c | ≤ max | max 0.25 (coarse limit; A&M printed 156 citing Ref. 2) |
| Khalil WDD-mimic | balance/cell-relation mimic | 0.2247c | monotone ↓ | 0.131 (A&M printed 156 citing Ref. 3) |
| **M4S (35)** | four-step with C lagged ≡ direct DFE of diffusion | 0.2247c (S16); 0.0 (S2) | **peak 0.67c (S2) → 0.50c (S16) at ~1–2 mfp** | **→ 0** |

Notes: (a) the 0.25 A&M quote and the memo §2.6 K_N ≤ 0.300 are BOTH Larsen-paper
numbers but different suprema (K_N table = sup over the α model family at each N;
0.25 = the coarse-mesh LD limit) — don't conflate. (b) A&M also analyzed the
**θ = 1 LD variant "in the notation of Ref. 15"** — Ref. 15 = **Larsen & Morel 1989
JCP 83:212, the ORPHEUS slab-LD reference** (issue #158 record) — and found
"performance essentially identical to that of the standard method" (printed 159,
verified). So the M4S numbers transfer to the LM-1989 LD family ORPHEUS implements.
(c) For the R5-gated ORPHEUS LD arm the trade is now explicit: fully-consistent
Larsen (78)/(81) buys the tighter uniform bound (≤ ~0.25–0.30c) at slab-only,
algebraically heavy cost; M4S buys always-eliminable simplicity, multi-D reach, and
the thick-cell ρ → 0, at the price of a ~2–3× higher intermediate-thickness peak and
the A&L-recorded divergence on skewed tetrahedral cells (memo §4.2 [300]) — M4S's
stability record is orthogonal-grid.

### 7.5 Stability results (§V–VI, printed 155–159, PDF 12–16)

**Fourier machinery** (verified PDF 12–13): subtract converged equations → error
equations (37) [LD transport: balance (37a), slope (37b), moments (37c), upstream
closure (37d)] + (38) [low-order] + (39) [two-moment update]; ansatz (41);
⟹ 2×2 eigenproblem **(42) ωv = [I − L]⁻¹[H − L]v**. Search over c ∈ (0,1),
λh ∈ (0,π); worst case always c → 1 (the memo-wide pattern).

**Slab LD numbers** (Fig. 2, printed 156, verified on the render):

- fine mesh: ρ → analytic DSA values — **0.0 for S2** (slab S2 ≡ P1: the K_2 = 0
  exactness analog survives the inconsistency), ≈ **0.2247 for S16**;
- intermediate (~1–2 mfp): peak **0.67c (S2), ~0.55c (S4/S8), 0.50c (S16)**;
- thick: **ρ → 0 monotonically** for all quadratures.

**2-D bilinear DFE** (§V.B + App. B; Tables II–IV, printed 157, PDF 14, verified via
text layer + render spot-check): 4×4 eigenproblem (four unknowns/zone), parameters
(quadrature, λ_x h_x, λ_y h_y, c, σ_t h_x, σ_t h_y); global-extremum search flagged
"fairly difficult … sometimes very sharp" peaks. S4 level-symmetric table (Table III,
ρ at c → 1; rows/cols = σ_t h_x / σ_t h_y ∈ {10⁻³, 10⁻², 0.1, 1, 2, 5, 10, 10², 10³},
lower triangle):

| σth_y \ σth_x | 0.001 | 0.01 | 0.1 | 1.0 | 2.0 | 5.0 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 0.25 | | | | | | | | |
| 0.01 | 0.25 | 0.25 | | | | | | | |
| 0.1 | 0.23 | 0.25 | 0.23 | | | | | | |
| 1.0 | 0.49 | 0.49 | 0.49 | 0.49 | | | | | |
| 2.0 | 0.55 | 0.55 | 0.55 | 0.55 | 0.55 | | | | |
| 5.0 | 0.44 | 0.44 | 0.44 | 0.49 | 0.55 | 0.44 | | | |
| 10 | 0.28 | 0.28 | 0.28 | 0.49 | 0.55 | 0.44 | 0.28 | | |
| 100 | 0.18 | 0.18 | 0.16 | 0.49 | 0.55 | 0.44 | 0.28 | 0.03 | |
| 1000 | 0.18 | 0.18 | 0.16 | 0.49 | 0.55 | 0.44 | 0.28 | 0.03 | 0.003 |

S2 (Table II) peaks 0.67 at (2,2) mfp; S8 (Table IV) peaks 0.50 at (2,2); both → ~0
in the both-thick corner; fine-limit analytic values 0.500 (S2) / 0.228 (S8) from
their Ref. 14 (M.L. Adams, "…The Fine-Mesh Limit," in preparation at the time).
Everything < 1 everywhere searched ⟹ the paper's unconditional-stability claim (2-D
orthogonal grids). ⚠ The high-aspect-ratio off-diagonal entries (e.g. 0.49–0.55
persisting at σ_t h_x = 1–2 for ANY σ_t h_y) foreshadow exactly the aspect-ratio
degradation A&L §IV.D tabulates for partially-consistent families.

**Numerical evidence** (§VI, printed 158–159, verified): estimator (44)/(45)
(successive-difference ρ, Δx-weighted L2 — same family as memo §6.3's A&L
(1.15)–(1.20) row).

- Test 1 (Fig. 3): homogeneous slab/sphere, σ_t = 1, σ_s = 0.99 (c = 0.99), S = 1
  left half / 0.01 right half, reflecting-left + vacuum-right, S8. Table V (verified
  on render): mesh 0.01 mfp → 7 iters to 10⁻⁶, ρ̂ ≈ 0.203/0.205 (slab/sphere);
  2.0 mfp → 13 iters, ρ̂ ≈ 0.471 (the intermediate peak); 10 mfp → 10 iters,
  ρ̂ ≈ 0.250; 100 mfp → 4 iters, ρ̂ ≈ 0.017. **Slab and sphere columns identical to
  ~±0.005 at every point.**
- Test 2 (Fig. 5): the "iron-water" problem — 4 regions (water/water/iron/water),
  σ_t = 3.3333/3.3333/1.3333/3.3333, σ_s = 3.3136/3.3136/1.1077/3.3136, source
  S = 1 in region 1 only, boundaries at 12/15/21/30 cm, reflecting-left +
  vacuum-right — the SAME problem family as the memo §3.1 McCoy-Larsen shielding
  spec (A&M cite it as "iron-water test problem²,³"). Table VI (S8 Gauss-Legendre,
  verified): ρ̂ from 0.19 (1 cell/region) peaking 0.46 (4 cells/region) back to 0.19
  (200 cells/region); every entry BELOW the model-problem prediction for the
  worst-case cell thickness; sphere ≡ slab again.

**The empirical upper-bound claim** (§VI intro + conclusions, printed 158/159): the
infinite-medium Cartesian Fourier predictions "appear to be upper bounds for slab and
spherical geometries, heterogeneous media, and variable mesh spacings."

**Sphere caveat** (printed 158): "we cannot Fourier-analyze the spherical-geometry
method" — the curvilinear arm is empirically gated only, exactly the memo §6.6.5
posture; A&M is the empirical precedent (sphere ρ̂ ≡ slab ρ̂ on every test).

### 7.6 Multidimensional + curvilinear statements

- **Any grid / any coordinate system** via the first derivation (§IV.A weak form
  (18)): the recipe is coordinate-free (surface integral of upstream partial currents
  + volume Galerkin terms).
- **Spherical** (Appendix A, printed 160–162, PDF 17–19): LD DFE spherical transport
  (A.1)–(A.4) — includes the angular-redistribution (γψ)_{m±½} terms with the ⅙/⅓
  slope-weighted redistribution moments and volume moments V_k, W_k, X_k — and the
  matching DFE diffusion (A.6) + (A.9) (update (A.8)), tridiagonal-formable
  ((A.11)–(A.14)). This is the corpus's ONLY worked curvilinear DFE-DSA low-order
  operator (Morel §5.5 gives the curvilinear D correction for cell-centered schemes,
  not DFE).
- **2-D Cartesian bilinear** (Appendix B, (B.1)–(B.10)): rectangular cells, four
  unknowns/zone — the operator behind Tables II–IV.
- **Eigenvalue + multigroup reach** (printed 160, verified): "our DFE diffusion
  discretization can be used to accelerate DFE transport eigenvalue iterations¹⁹ and
  group-to-group scattering iterations as well as within-group iterations, in the
  same way that any linear DSA scheme can." Ref. 19 = Adams & Larsen, "Linear
  Diffusion-Synthetic Acceleration of k-Eigenvalue Problems," Proc. Santa Fe 1989 —
  a named LINEAR k-DSA primary for §6.6.3 (previously the memo carried only the
  Gelbard-et-al. [81] thread from A&L §VIII.E).

### 7.7 Notation map (A&M → ORPHEUS/memo)

- **Σw = 1** (proof: (25a) half-range Σw(n·Ω)² = 1/6, (25b) α ≈ ¼) — the §6.5
  Alcouffe/Larsen/McCoy/Morel camp. Map before comparing with A&L (Σw = 2) rows.
- Their **α** = discrete half-range current weight (≈ ¼) = Larsen's discrete γ_N
  (§2.4). ⚠ COLLISION: NOT Larsen's WD spatial α (§2.2), NOT Morel's albedo α (§5.4),
  NOT the BMC curvilinear α. Four α's now live in this memo — qualify every use.
- Their **β**: 0 = reflective, 1 = vacuum/incident (γ ≡ 1 − β in the f-form) —
  reverse direction of Morel's albedo.
- Their **C(r_k)** = the first-moment interface-jump surface term — the M4S-specific
  object; no analog in §§1–5 (fully-consistent reductions keep it live).
- Iterates l+½ / l+1 = the memo convention exactly; their f₀ = the memo's f
  (correction form); LD pair (f₀k, f₀k^x) = cell-average + slope moment = the
  ORPHEUS LD moment pair (Larsen-Morel 1989 family, θ = 1 variant covered).
- Their "standard linear DSA" = the memo §2 four-step (their Ref. 2 = Larsen Part I).

### 7.8 What M4S changes in the memo's picture

1. **"Consistency is sufficient but NOT necessary" now has its existence proof
   locally** (their conclusions, printed 159–160: "strict consistency … is not
   necessary. However, we still do not know how much inconsistency a scheme can
   tolerate"). Read together with McCoy-Larsen §3.3 (Table II divergence of the
   WD-transport/DD-low-order mix): *partial consistency is not one bucket* — WHICH
   terms are lagged/mismatched decides stability. Lagging the first-moment interface
   jumps C (M4S) is benign on orthogonal grids; grafting a foreign low-order operator
   + BCs (Table II) is fatal. The memo §6.1's "derived-by-moment-reduction is the
   operative property" survives intact — M4S IS a moment reduction, just with one
   extra lagged block.
2. **The §4.2 taxonomy row for M4S is now firsthand-verified** (was secondhand via
   A&L [169]): rapid convergence + unconditional stability demonstrated on orthogonal
   grids (1-D/2-D Fourier + slab/sphere numerics); the skewed-tet divergence remains
   A&L/Warsa-sourced (not in this paper — its record is orthogonal-grid only).
3. **Thick-cell ρ → 0** is a property NO fully-consistent scheme in the memo has
   (Larsen DD/LD saturate at ~0.25c–0.3c): relevant if ORPHEUS ever needs
   radiative-transfer-style thick cells on the LD arm.
4. **§6.6.8 (LD-arm reference) gains a second published route** — see the §7.4 table
   and trade-off note.

---

## 8. Warsa-Wareing-Morel-McGhee-Lehoucq 2004 — Krylov (IRAM) k-eigenvalue iterations

J.S. Warsa, T.A. Wareing, J.E. Morel, J.M. McGhee & R.B. Lehoucq, "Krylov Subspace
Iterations for Deterministic k-Eigenvalue Calculations," NSE **147**(1), 26–42 (May
2004), DOI 10.13182/NSE04-1. LOCAL. Page map: printed ≈ PDF + 24 (printed 26 = PDF 2).

> ⚠ **Identity disambiguation (load-bearing)**: this is NOT the "degraded-DSA" paper
> the roadmap wanted. TWO Warsa-et-al. papers live in NSE 147 (2004): **147(1), 26–42**
> = THIS one (k-eigenvalue IRAM; five authors, + McGhee & Lehoucq) and **147(2),
> 218–248** = "Krylov Iterative Methods and the Degraded Effectiveness of Diffusion
> Synthetic Acceleration for Multidimensional S_N Calculations in Problems with
> Material Discontinuities" (three authors; STILL ABSENT locally). The sibling
> relationship is visible inside this paper: its Ref. [25] cites the absentee in
> submitted-manuscript form under the early title "Diffusion Synthetic
> Acceleration—Part I: Deficiencies in Multi-Dimensional Heterogeneous Problems"
> (verified, reference list printed 41). Do not let a citation to "Warsa et al. 2004"
> resolve to the wrong member.

> **Provenance note (rev. 2)**: same Mistral-key 401 as §7 — search surface is the
> pypdf text-layer file `scratch/literature_ocr/<stem>.textlayer.md`; equations and
> ALL iteration-count/tolerance tables verified on rendered PDF pp. 4, 8–9, 11–12, 14
> (`Read pages=`). Regenerate the true sidecar when the key is fixed.

### 8.1 The operator the Arnoldi iteration runs over (§II.A–II.B, II.D; printed 27–29, PDF 3–5)

**Transport discretization** (Eqs. (3a)/(3b), printed 27–28, verified PDF 4):
multigroup 3-D SN, DFEM with linear basis on tetrahedra (4 unknowns/cell), upstream
interface fluxes (3b), **fully lumped including the surface integral** — the lumping
"preserves the diffusion limit in thick, diffusive regimes" (their Ref. 9 = Adams
2001 NSE 137:298 — the same thick-limit doctrine as memo §4.6). Vacuum (Γ = 0) and
specular-reflective BCs (reflective faces axis-aligned; quadrature contains the
reflected pairs Ω̂_m′ = Ω̂_m − 2n̂(Ω̂_m·n̂)).

**Operator form** (verified PDF 4): Eq. (4) **Lψ = MSDψ + (1/k)FDψ** — L transport,
S scattering, F fission, D angular→moments (φ = Dψ, "all the scalar flux moments
needed to construct the scattering source"), M moments→angular source. Rearranged to
the power iteration (5):

φ^{(ℓ+1)} = (1/k^{(ℓ)}) A φ^{(ℓ)},  **A = D H⁻¹ F,  H = L − MSD**

with k updated by the norm ratio (6) (fission-rate norm (7) ‖φ‖_F = Σ_g Σ_s νσ_f
φ⁰₀); the Rayleigh-quotient alternative (8)/(9) "can usually improve the efficiency"
of PI. The k-eigenproblem handed to the eigensolver is the STANDARD problem (10)
**Aφ = kφ on the moment vector** — the transport inverse H⁻¹ lives INSIDE the
operator whose action each Arnoldi step computes. BCs "ignored without loss of
generality to facilitate discussion" (so H carries no explicit −B term in print).
The 1/k^{(ℓ)} scaling in (5) is deliberate: prevents over/underflow AND lets PI
warm-start inners from the previous iterate (§II.E, printed 31).

**The eigensolver** (§II.D, printed 29–30): IRAM as implemented in ARPACK (their
Refs. [2] Sorensen 1992 SIAM J. Matrix Anal. Appl. 13:357; [7] Lehoucq-Sorensen-Yang
ARPACK Users' Guide 1998). Mechanics recorded: Krylov subspace K_m(A, φ₀); Ritz
pairs by the Galerkin condition; m×m upper-Hessenberg projection + QR; Ritz
(backward-error) estimates for free; implicit restart = implicitly-shifted QR with
ARPACK's exact shifts (the m − r smallest Hessenberg eigenvalues), "mathematically
equivalent to computing a new Krylov space with a starting vector that is a linear
combination of the r dominant Ritz vectors"; convergence rate (11) = the
polynomial-filtered ratio (power iteration recovered at r = 1, P(x) = x; subspace
iteration at all-zero shifts); ≡ accelerated subspace iteration [23]. Implementation
via reverse communication, "wrapped around the existing power iteration" — the
abstract's "only modest changes to existing power iteration coding" claim.

**Stopping** (§III preamble, printed 32, verified PDF 8): ARPACK backward error (14)
‖Ax − θx‖₂ ≤ |θ|ε via Ritz estimates; the PI-side identity (15) ‖r^{(ℓ)}‖₂ =
|k^{(ℓ)}|·‖φ^{(ℓ+1)} − φ^{(ℓ)}‖₂ justifies the matched PI test ‖φ^{(ℓ+1)} −
φ^{(ℓ)}‖₂ ≤ ε. ε = 10⁻⁴ (outer) for both methods; convergence is judged on the
EIGENVECTOR. The "k only" rows of Tables IV/VI/VII (eigenvalue to 10⁻⁵ absolute)
show the PI eigenvalue converging ~2× sooner than its eigenvector (symmetric-case
theory: quadratic vs linear in δ [36]) — "a property the IRAM does not share"
(printed 35).

**Spectrum caveats** (printed 28–29): k real proven only for MONOenergetic SN
(their [12] Modak-Sahni-Paranjape ANE 22:359); multigroup unproven but empirically
the 10–20 largest eigenvalues are real across their problems; the DFEM-on-tets
discretization DESTROYS the symmetrizability of Refs. [13–15] (Ramone-Adams-Nowak;
Sanchez-Santandrea). Power iteration is compared UNACCELERATED because no robust
Chebyshev parameter set could be found and simple SOR diverged (§II.C printed 29 +
§III.B.1 printed 36) — nonsymmetry makes classical spectrum-guess acceleration
fragile, which is itself part of the paper's robustness argument.

### 8.2 The inner-solve architecture and BOTH DSA postures (§II.E + §III preamble; Tables I/II, printed 31–33, verified PDF 8–9)

**Energy nesting** (§II.E, Eqs. (12)/(13), printed 30–31): split S = S_L
(downscatter + within-group, block lower triangular) + S_U (upscatter). No upscatter
⟹ the action of H⁻¹ = block forward substitution down the group ladder, one
within-group transport solve per group. With upscatter ⟹ an INTERMEDIATE iteration
level (block Gauss-Seidel — Eq. (13) is "just an outer Gauss-Seidel iteration" — or
FGMRES(20) over the coupled thermal block; Fig. 7's 2×2 partition: forward-substitute
the downscatter-only block, iterate the 4-group thermal block). Nesting under IRAM:
**IRAM outermost → intermediate H⁻¹ (upscatter) → within-group inners**. Under PI the
H⁻¹ Gauss-Seidel is what tradition calls "the outers."

**The two DSA postures, side by side** (printed 32, verified verbatim + Tables I/II):

- **Power-iteration arm = SI+DSA inners**: "The power iteration results are computed
  using DSA-accelerated source iteration, and solutions from the most recent outer
  iteration are used as initial guesses" — Table I labels the level "DSA (inner
  acceleration)".
- **IRAM arm = DSA-preconditioned-Krylov inners**: "The ARPACK results are computed
  using a preconditioned Krylov method for the inner iterations, using the source
  vector (the uncollided flux) as the initial guess [32]" — Table II labels the level
  "DSA (inner preconditioner)". Krylov flavor: FGMRES(10) (cylinder problem) /
  BiCGStab (MOX problem); flexible GMRES because the DSA-CG preconditioner is itself
  an inner iteration (variable preconditioner).

**Which DSA**: in BOTH postures it is "the **partially consistent, simplified DSA
method of Wareing, Larsen, and Adams [33]**" (1991 Pittsburgh proceedings): solve a
linear CONTINUOUS-finite-element diffusion equation by conjugate gradients, then map
back to update the discontinuous fluxes [27, 33]. It is NOT the fully consistent
DFEM DSA the same authors built in [27] (= Warsa-Wareing-Morel, "Fully Consistent
Diffusion Synthetic Acceleration of Linear Discontinuous S_N Transport
Discretizations on Unstructured Tetrahedral Meshes," NSE **141**, 236 (2002),
verified in the reference list) — the production choice here is the cheap
partially-consistent operator, with the degraded-effectiveness analysis delegated to
the absent sibling [25]. Note the vocabulary: by 2004 "partially consistent" is the
standard label (it is NOT in Adams-Martin 1992 itself — see §7.2).

**Tolerance cascade** (Tables I/II, verified on the render — the text layer garbles
the nesting):

- PI (Table I): outer 10⁻⁴; inner min(10⁻², 10‖r_outer‖), floor 10⁻⁵ when
  ‖r_outer‖ = 0; DSA 10⁻¹‖r_inner‖, floor 10⁻⁶.
- IRAM (Table II): outer 10⁻⁴; each lower level's tolerance is
  max(10⁻¹‖r‖, 10⁻¹ min(10⁻¹, 10⁻¹‖r‖/min(1, ‖r‖))) in the residual ‖r‖ of the level
  above (r = r_outer for upscatter, r_o for inners, r_inner for DSA; r_o ≡ r_upscat
  when upscatter is present else r_outer); floors when ‖r‖ = 0: upscatter 10⁻⁵,
  inner 10⁻⁵/10⁻⁶ (without/with upscatter), DSA 10⁻⁶/10⁻⁷.

This is the **BFG (Bouras-Frayssé-Giraud) inexact-matvec relaxation** [28–30], with
Simoncini-Szyld theory [31] fixing the inverse-proportionality constant (they use
0.1): tolerances START STRICT while the first Krylov subspace of dimension ncv is
built, then RELAX as the outer residual falls. Architecturally this REPLACES
warm-starting: Arnoldi basis vectors are not flux-like, so "good initial guesses are
not available for subsequent inner iterations in the IRAM iterations" (printed 31) —
the mirror image of PI's warm-start economics (Figs. 5 vs 6: PI inners decay to 1–2
per outer late in the iteration; IRAM inners start expensive during subspace
construction, then decay as tolerances relax).

**Initialization**: random start; optionally init = 5 power iterations, each with a
SINGLE inner per group — cheap, and "certainly worth the cost" (printed 33–34).

### 8.3 Robustness evidence vs dominance ratio (§III; Tables IV–VII, verified PDF 8–9/11–12/14)

**Cylinder family** (§III.A, printed 32–33): r = 3.5 cm × 9 cm cylinder; 5-cm B-10
absorber sandwiched between 1-cm water and 1-cm HEU layers; vacuum everywhere;
13,500 tets; five 0.1-cm disks beside the absorber flipped water→boron one at a time
to break the mirror symmetry: δ = |k₂|/|k₁| ∈ {0.970, 0.977, 0.984, 0.990, 0.996,
0.999(symmetric)}; k barely moves, the eigenvector/δ move. 5-group Hansen-Roach
collapse (16→5, map printed 33); HEU ρ=19.0 (95/5 U-235/U-238), water, B-10 ρ=10.0;
S4 triangular Chebyshev-Legendre; FGMRES(10) inners. Fig. 2: the two largest-δ PI
runs fail to reach 10⁻⁴ within the 2500-iteration plot window.

Table IV — outer iterations (verified; (nev, ncv, init) rows):

| method | δ=0.970 | 0.977 | 0.984 | 0.990 | 0.996 | 0.999 |
|---|---|---|---|---|---|---|
| IRAM (1,10,0) | 20 | 20 | 20 | 20 | 20 | 20 |
| IRAM (1,10,5) | 15 | 15 | 15 | 15 | 15 | 15 |
| IRAM (1,5,0) | 20 | 20 | 20 | 20 | 23 | 23 |
| IRAM (1,5,5) | 14 | 14 | 14 | 14 | 14 | 17 |
| IRAM (2,10,0) | 18 | 18 | 18 | 25 | 26 | 26 |
| IRAM (2,10,5) | 18 | 18 | 18 | 18 | 18 | 18 |
| IRAM (2,5,0) | 20 | 20 | 20 | 20 | 23 | 25 |
| IRAM (2,5,5) | 14 | 14 | 14 | 14 | 14 | 17 |
| Power iteration | 468 | 602 | 833 | 1329 | 3170 | 8824 |
| PI, k-only (10⁻⁵) | 229 | 287 | 388 | 586 | 1206 | 1739 |

Table V — total inners (verified; sample rows): IRAM (1,5,5): 327/327/327/332/322/372;
IRAM (1,10,0): 460/460/458/457/447/447; worst IRAM row peaks 550; **PI:
2503/3124/4218/6691/15892/44156**. Inners are 95–99% of FLOP cost for both methods,
so inner counts ≈ total cost; Fig. 3: IRAM is 1–2 orders of magnitude cheaper than
unaccelerated PI, insensitive to nev/ncv. **The headline robustness fact: IRAM outer
counts are FLAT in δ (14–26 across δ = 0.970→0.999) while PI grows ~19× — and the
second eigenpair (nev = 2) costs almost nothing extra.**

**MOX C5G7 benchmark** (§III.B, printed 35–38): the OECD/NEA C5G7 MOX 3-D spec
(their [37], Lewis NEA/NSC/DOC(2001)17 App. C); 954,527 tets / 169,745 vertices /
~28M DOF (serial!); 7-group transport-corrected isotropic; S4 triangular C-L;
BiCGStab inners; reflective on the three coordinate planes, vacuum elsewhere.

- Without upscatter (III.B.1, upscatter entries zeroed; δ ≈ 0.985), Table VI
  (verified): IRAM 25–44 outers / 1050–1578 inners (best (1,10,5): 25/1050;
  (1,5,5): 32/1132) vs PI 860/6144 (k-only 296/2196) ⟹ **IRAM ≈ 5× cheaper**; wall
  clock ~23–34 h vs ~150 h.
- With upscatter (III.B.2, δ ≈ 0.987), Table VII (verified): IRAM (1,5) 32 outers /
  6011 inners / 504 Gauss-Seidel upscatter iters (FGMRES(20) upscatter variant:
  32/6991/309; (2,10): 33/6030/483) vs PI 1144/8144 ⟹ IRAM ≈ **75% of PI cost**
  only; wall ~143–157 h vs ~203 h. The nested upscatter level eats the advantage —
  the conclusions' declared "significant drawback … the inability to efficiently
  treat problems with energy upscatter" (printed 39, PDF 16).

Future-work list (conclusions, printed 38–39): in-line IRAM (reuse transport data
structures), combining Chebyshev/coarse-mesh-rebalance, shifted-inverse iteration,
parallel P_ARPACK [38].

### 8.4 What this adds to the memo (consumer: §6.6.3 k-eigenvalue coupling)

1. **A linear, parameter-free k-loop that leaves DSA untouched**: the dominance-ratio
   robustness problem is solved at the OUTER level (Arnoldi subspace, exact shifts,
   no user-estimated spectrum parameters), fully decoupled from HOW each transport
   application is solved. DSA remains a within-group device in EITHER posture
   (§8.2). This is the missing third option for §6.6.3 next to (a) Alcouffe's
   nonlinear D̂-correction lineage and (b) Gelbard-style linear DSA-eigenvalue
   variants: (c) keep the k-loop linear-algebraic (power → Arnoldi on A = DH⁻¹F)
   and let DSA precondition the inners.
2. **The inner-tolerance design fact**: an Arnoldi outer forfeits warm starts;
   the BFG relaxed-tolerance cascade (constant 0.1, strict-during-subspace-build)
   is the published replacement, with Simoncini-Szyld backing. Any ORPHEUS
   IRAM-style k-loop inherits this trade.
3. **Posture pairing precedent**: PI ↔ SI+DSA(acceleration)+warm-start; Arnoldi ↔
   DSA-preconditioned-Krylov(FGMRES/BiCGStab)+relaxed tolerances — the same pairing
   logic as memo §4.4's "GMRES-outside-DSA" quote, now demonstrated at production
   scale (28M DOF).
4. **What it does NOT contain**: no Fourier/spectral analysis of DSA, no consistency
   derivations, no degraded-DSA evidence (that is the ABSENT sibling NSE 147(2):218,
   here cited as [25] in submitted form). Its DSA content is purely architectural
   (which posture, which tolerances, WLA partially-consistent operator [33]).

### 8.5 Notation map (Warsa et al. → ORPHEUS/memo)

- Their **L** = streaming + collision + interface coupling (the sweepable block) ⟹
  ORPHEUS L + C (BCs suppressed in their operator discussion); their **H = L − MSD**
  = the within-group-coupled transport operator ⟹ the ORPHEUS honest A = L + C − S
  (mod B). Their **A = DH⁻¹F** = the moment-space projection of the ORPHEUS
  eigenvalue map K = A⁻¹F (memo §6.5 / #226 vocabulary): eigenproblem posed on
  MOMENT vectors, not angular flux — Arnoldi storage is per-moment, not per-ordinate.
- Their **D** (ψ → all scattering moments) and **M** (moments → angular source) are
  the moment-projection / harmonic-lift pair of the SN scattering machinery.
- Their fission term ((3a), verified): (1/k)χ_{g,s} Σ_{g′} νσ_{f,g′,s} φ⁰₀ — χ
  strictly LOCAL to the emission cell, the same indexing law as the memo's sources.
- **δ** = |k₂|/|k₁| (the text layer renders δ as "d" — OCR artifact); their ℓ = the
  memo's outer index l.
- "Outer" is REBOUND under IRAM: one outer = one application of A (one full H⁻¹
  sweep-solve down the group ladder), while under PI "outer" is the classic
  power/Gauss-Seidel step. Table headers use the IRAM sense.
- Quadrature: triangular Chebyshev-Legendre S4; cgs units; norm (7) = fission-rate
  norm; eigenvector-based stopping at 10⁻⁴ both methods (§8.1).
