# DSA literature memo — consistent Diffusion Synthetic Acceleration for the ORPHEUS SN solver (issue #2, Phase 3)

**Purpose**: the literature record a future implementing session needs to build
consistent DSA for the slab DD SN solver and gate it against published spectral-radius
predictions. Every claim below carries source + equation number + printed page (and PDF
page for `Read pages=`). Extracted 2026-07-26 by literature-researcher; all equations
spot-verified against the rendered PDF pages (the scan is the SSOT).

---

## 0. Source inventory

All five primary sources are LOCAL in `scratch/literature/`, each with a Mistral-OCR
sidecar in `scratch/literature_ocr/<stem>.md` (`## p. N` = 1-based PDF page).

| # | Paper | Journal | DOI | Local status | Page map (printed → PDF) |
|---|---|---|---|---|---|
| 1 | Alcouffe (1977), "Diffusion Synthetic Acceleration Methods for the Diamond-Differenced Discrete-Ordinates Equations" | NSE **64**(2), 344–355 | 10.13182/NSE77-1 | PDF + sidecar ✓ | printed ≈ PDF + 342 (printed 344 = PDF 2) |
| 2 | Larsen (1982), "Unconditionally Stable DSA Methods for the Slab Geometry Discrete Ordinates Equations. Part I: Theory" | NSE **82**(1), 47–63 | 10.13182/NSE82-1 | PDF + sidecar ✓ | printed ≈ PDF + 45 (printed 47 = PDF 2) |
| 3 | McCoy & Larsen (1982), "… Part II: Numerical Results" | NSE **82**(1), 64–70 | 10.13182/NSE82-2 | PDF + sidecar ✓ | printed ≈ PDF + 62 (printed 64 = PDF 2) |
| 4 | Adams & Larsen (2002), "Fast Iterative Methods for Discrete-Ordinates Particle Transport Calculations" | PNE **40**(1), 3–159 | 10.1016/S0149-1970(01)00023-3 | PDF + sidecar ✓ + **FULL prior extraction** `.claude/plans/phase_i_survey_adams_larsen_2002.md` | journal p. N = PDF N − 2 |
| 5 | Morel (1982), "A Synthetic Acceleration Method for Discrete Ordinates Calculations with Highly Anisotropic Scattering" | NSE **82**(1), 34–46 | 10.13182/NSE82-A19026 | PDF + sidecar ✓ | printed ≈ PDF + 32 (printed 34 = PDF 2) |

Note that #2, #3, and #5 all sit in the same NSE issue (82(1), 1982) — the DSA issue.

**ABSENT from the local folder** (roadmap named them; NOT acquired online per the
delegation rule — user to add if wanted):

- **Adams & Martin (1992)**, "Diffusion Synthetic Acceleration of Discontinuous Finite
  Element Transport Iterations," NSE **111**, 145 — the M4S paper. Coverage here is
  secondhand only, via A&L 2002 [169] (§IV.D: M4S is rapidly convergent on orthogonal
  grids, DIVERGES on skewed tetrahedral cells per Warsa et al. [300]).
- **Warsa, Wareing & Morel (2004)** (Krylov-DSA robustness; the "fully consistent DSA
  can degrade in multi-D heterogeneous problems, Krylov fixes it" paper). A&L 2002
  predates it; the review carries only the precursor thread (consistent-DSA
  heterogeneity degradation to σ ≈ 0.88 with "no known remedy", pp. 69/139, and the
  2001 M&C proceedings version [300]).

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
A **complete** prior extraction exists: `.claude/plans/phase_i_survey_adams_larsen_2002.md`
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
3. **Eigenvalue coupling**: Part I/II are fixed-source; Alcouffe's k-eigenvalue
   variants are nonlinear (D̂/removal correction, DANTSYS lineage) with fixup
   interactions; the linear alternative (Gelbard et al. [81], A&L §VIII.E) costs an
   adjoint low-order eigenproblem. ORPHEUS #2 scope = within-group/fixed-source
   acceleration first; the k-loop design is a separate decision.
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
7. **Absent primaries**: Adams & Martin 1992 (M4S; needed for a DG/LD multi-D arm);
   Warsa-Wareing-Morel 2004 (Krylov robustness for degraded consistent DSA);
   Yavuz-Larsen 1988 (reflective BCs). User to acquire if/when those arms open.
8. **LD-arm reference** (gated): Larsen Part I §V — LD α_m = ε_m/(3+|ε_m|),
   tridiagonal (78) with (81) coefficients, K_N ≤ 0.300; plus the A&L p. 73
   presentation-dependence warning.

### 6.7 Reading order for the implementing session

A&L §§II.B + III.B (tutorial + consistency mechanics, with the survey file) →
Larsen Part I §§II–III (the recipe + BCs — the 3a reference) → McCoy-Larsen (the gate
tables) → Alcouffe (DD specifics, 2-D corner scheme, eigenvalue variants; mind §1.5
sign errata) → Morel (P1 extension + curvilinear + albedo rows) — matching the survey
Task-5 recommendation.
