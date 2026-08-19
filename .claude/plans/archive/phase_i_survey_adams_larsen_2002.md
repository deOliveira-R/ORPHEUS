# Phase-I survey — Adams & Larsen (2002), "Fast Iterative Methods for Discrete-Ordinates Particle Transport Calculations"

**Full citation**: Marvin L. Adams and Edward W. Larsen, "Fast Iterative Methods for
Discrete-Ordinates Particle Transport Calculations," *Progress in Nuclear Energy* **40**(1),
pp. 3–159 (2002). PII `S0149-1970(01)00023-3` (printed on p. 3; DOI
`10.1016/S0149-1970(01)00023-3`). Provenance is direct: all claims below are extracted from
the local artifact, not from secondary citations.

**Source artifact**: `scratch/literature/Adams-Larsen(2002) Fast iterative methods for
discrete-ordinates particle transport calculations.pdf` (OCR'd; 158 PDF pages). Page
mapping: journal p. *N* = PDF page *N*−2 (1-based). All page numbers below are **journal
pages**. Extraction notes: OCR text is reliable for prose; every load-bearing equation below
was transcribed from the **rendered PDF page**, not the OCR (the OCR renders ⅓Λ² as "½A2"
in Eq. (3.41)/(3.43) — trust the transcriptions here, which come from the page images).

**Survey context**: fills corpus-plan gaps **M4b** (convergence theory — textbooks have
zero: no `spectral radius|Fourier analysis|convergence rate` hits, no ρ ≈ c predictor) and
**M2** (symptom→cause diagnostics). This review is where that missing content lives — it is
THE iteration/acceleration reference for the S_N docs' planned convergence/admissibility
section.

---

## Task 1 — §3.2 table row

| Column | Verdict | Evidence |
|---|---|---|
| **S_N page count** | **157 pp, all S_N** (pp. 3–141 content; 142–158 references ≈ 300 entries; 158–159 acronym appendix) | Scope statement, p. 6: the review covers iterative methods for the discrete-ordinates discretization only — "We do not discuss iterative methods for other angular discretizations, such as the spherical harmonics or P_N approximation, which possess different mathematical properties and difficulties." |
| **Cylinder?** | **≈ ⅓ page, zero equations, zero theory** (§IV.C, p. 80). The entire 1-D cylindrical treatment is two paragraphs: geometry definition + "apply slab lessons, it works." Contrast: 1-D **spherical** gets 3 pages (§IV.B, pp. 77–80) including a *new* continuous-geometry Fourier result (see Task 2.9). For discretized curvilinear S_N, no Fourier analysis exists at all (p. 79). | p. 80, verbatim: "As in spherical geometry, the practical development of acceleration schemes in cylindrical geometry has proceeded by applying the lessons learned from successful Cartesian approaches. As in spherical geometry, this strategy has been effective; cylindrical geometry has not yielded any major surprises [144, 175]." ([144] = ONEDANT manual; [175] = Palmer's 1992 Michigan PhD thesis.) |
| **Adjoint of S_N?** | **Absent as a topic.** 8 "adjoint" hits total, all instrumental: SAAF equation name (p. 47); self-adjointness of the CG operator/inner product (p. 109); a zero-eigenvalue existence argument via the adjoint of the upscattering matrix (p. 115); Gelbard et al.'s linear k-eigenvalue scheme "requires the solution of one adjoint low-order eigenvalue problem" (p. 137); ref [277] title (Anistratov, adjoint QD, 2000); acronym list (p. 159). No adjoint transport equation is ever posed or iterated. | grep count over full text: `adjoint` → 8 hits at pp. 47, 109 (×2), 115, 137, 156 (×2 ref titles), 159. |
| **Eigenvalue posed independently?** | **YES — own chapter** (Ch. VIII, pp. 126–137, 12 pp), with its own spectral theory: k-eigenvalues are eigenvalues of A = L⁻¹P (Eq. (8.5)); the **dominance ratio** r ≡ \|k₂\|/k₁ (Eq. (8.6)) is the eigenvalue-iteration analog of c; nested outer/inner structure explicit; **shifted power iteration** (§VIII.B, Eqs. (8.18)–(8.25)) with r_SPI = \|κ−k₁\|/\|κ−k₂\| · \|k₂\|/\|k₁\|; Chebyshev-for-k (§VIII.C); nonlinear QD-for-k (§VIII.D); Alcouffe's nonlinear k-DSA + a linear alternative (§VIII.E). | p. 127, verbatim: "The dominance ratio r plays the same role for eigenvalue problems as the scattering ratio c does for fixed-source problems: if r is small, the simplest iteration scheme to determine k₁ (power iteration) converges rapidly, while if r is 'large' (near unity), power iteration converges slowly." Generality claim p. 130: "The above results for the PI and SPI schemes are very general and hold for transport problems in general geometries, with anisotropic scattering, and with (multigroup) energy dependence." |
| **Verification?** | **`benchmark` → 1 hit** (p. 95: a 3-D IAEA reactor benchmark used to test ICSA); **`manufactured` → 0; `verification` → 0.** The review's de-facto verification discipline is instead pervasive and specific: **Fourier-predicted spectral radius vs. observed iteration behavior** (`Fourier` → 66 hits, `spectral radius` → 69 hits), plus iteration-count tables (Table 1 p. 110, Table 2 p. 113, Alcouffe's three test problems p. 136). §III.J (pp. 68–69) is an explicit essay on the predictivity of the model-problem Fourier analysis for realistic problems. | p. 90: the linearized-rebalance Fourier predictions "agree very closely with numerical experiments of both the linearized and nonlinear methods [127, 150]"; p. 76: "the spectral radii approach the expected analytic value (0.2247c) for fine meshes"; p. 117 (upscattering): "essentially identical to the behavior predicted by the Fourier Analysis"; p. 68 (III.J): "the power of the Fourier analysis is not so much that it accurately predicts the behavior of transport iteration schemes in idealized problems, but that it also predicts the behavior in realistic non-model problems as well." |

---

## Task 2 — M4b DEEP harvest: the convergence theory

### 2.1 The Fourier-analysis METHOD (the teachable recipe)

The review's principal analysis tool, introduced §II.A (pp. 24–26) and restated for every
scheme. The recipe, exactly as the review executes it:

1. **Model problem**: infinite homogeneous medium (constant Σ_t, Σ_s), no boundaries;
   for discrete analysis add: uniform mesh h_j = h, solution bounded as j → ±∞ (p. 50).
2. **Error equation**: subtract the iteration equation from the exact equation; the equation
   for the error f^(ℓ) = ψ − ψ^(ℓ) is the iteration equation with Q = 0 (Eqs. (2.6)–(2.7),
   p. 24). Equivalently "set the source Q to zero and interpret the unknown functions as
   iteration errors" (p. 78).
3. **Fourier ansatz**: expand the error in modes. Continuous (Eqs. (2.39)–(2.42), p. 28–29):

   φ^(ℓ)(x) = ω^ℓ(λ) e^{iΣ_t λx},  ψ^(ℓ+1/2)(x,μ) = ω^ℓ(λ) α(λ,μ) e^{iΣ_t λx},
   φ^(ℓ+1/2)(x) = ω^ℓ(λ) β(λ) e^{iΣ_t λx},  F^(ℓ+1)(x) = ω^ℓ(λ) γ(λ) e^{iΣ_t λx}.

   NOTE the scaling: the wave number λ is **dimensionless, in mean-free-path units**
   (the phase is Σ_t λ x). Discrete (Eqs. (3.18)–(3.22), p. 50): same, evaluated at x_j and
   x_{j+1/2}, with per-direction edge/average amplitudes a_n, b_n.
4. **Reduce**: linear independence of modes ⇒ a per-λ algebraic system for the amplitudes;
   solve for the **iteration eigenvalue** ω(λ) (the per-iteration attenuation factor of mode λ).
5. **Spectral radius**: σ = sup_λ |ω(λ)|. Convergence iff σ < 1; the asymptotic error
   reduction per iteration is σ.

Foundational facts the analysis establishes (p. 26, four numbered observations):
modes are **independent** ("the Fourier modes for different λ are independent, and to
determine ω(λ), it is only necessary to consider the single Fourier mode"); the rate for a
general initial error is sup|ω|; and — the physical heart of the whole subject —

> "the maximum value of ω occurs for λ ≈ 0. The λ ≈ 0 modes correspond to long
> wavelengths; these are error modes that span large optical distances and have weak
> spatial gradients. By Eq. (2.20), the angular dependence of these modes is also weak
> … Thus, in the SI scheme, the most slowly converging component of the iterative
> solution has weak spatial and angular dependence. Source Iteration efficiently
> suppresses the error modes with strong spatial and angular variation (|λ| ≫ 0), but not
> the error modes with weak spatial and angular variation (λ ≈ 0)." (p. 26; Eq. (2.21):
> 1/(1+iλμ) = 1 − iλμ + O(λ²).)

This is the design principle for every accelerator in the review: the low-order operator
must be accurate exactly where the sweep is weak — nearly-flat-in-angle,
nearly-flat-in-space error modes (p. 27).

**Curvilinear caveat** (load-bearing for ORPHEUS): the separable ansatz fails in curved
geometry — see 2.9.

### 2.2 SI spectral radius: ρ = c — exact statements and assumptions

**Continuous transport, slab** (§II.A). Error attenuation (Eqs. (2.10)–(2.12), p. 25):

  A^(ℓ+1)(λ,μ) = (c/2) · 1/(1+iλμ) · ∫₋₁¹ A^(ℓ)(λ,μ′)dμ′,  c ≡ Σ_s/Σ_t,

  ω(λ) = c ∫₀¹ dμ/(1+λ²μ²) = (c/λ) tan⁻¹λ.   (2.12)

  σ_SI = max_λ ω(λ) = ω(0) = c.   (2.17, p. 26)

> "Thus, the error in the ℓ-th iterate decreases as c^ℓ. If c is small, the convergence is
> rapid, but as c approaches unity, the convergence becomes arbitrarily slow." (p. 26)

Assumptions: infinite homogeneous medium, isotropic scattering, one group, c < 1,
bounded iterates. Physical reading (p. 24): the convergence rate is "directly linked to the
mean number of collisions in a neutron lifetime"; ψ^(ℓ) with zero initial guess is exactly
the flux of particles that have scattered at most ℓ−1 times (Eq. (1.10), p. 8).

**Discretized (DD S_N), slab** (§III.A). With the discrete ansatz, Eqs. (3.23)–(3.24)
(p. 50) give a_n = (c/2)/(1 + iμ_nΛ) (3.27), where

  Λ ≡ (2/(Σ_t h)) tan(Σ_t h λ/2),   (3.28, p. 51)

  ω_SI(Λ) = (c/2) Σ_{n=1}^{N} w_n/(1+μ_n²Λ²),   (3.29)

  σ_SI = sup_Λ |ω_SI(Λ)| = ω_SI(0) = c.   (3.30)

> "This result holds for any even quadrature set order (N) and any spatial cell-size (h)."
> (p. 51)

I.e., **discretization does not change ρ_SI = c** — the discrete map Λ(λ) merely
reparametrizes the mode axis; Λ = 0 (flat mode) still attains c.

**Departures from the model problem** (§III.A, p. 51 — all three are docs-worthy
predictors):

- **Leakage**: in finite media "the observed spectral radius is somewhat less than c …
  The inclusion of leakage shortens neutron lifetimes and hastens the convergence …
  this effect diminishes as the physical system becomes increasingly optically thick."
- **Nonuniform grids**: "almost no change in the observed SI spectral radius, and in all
  cases the spectral radius is bounded from above by c."
- **Heterogeneous media** (spatially periodic analysis, ref [257] Zika & Larsen 1998; observed
  generally): for optically thick systems,

    min_x c(x) ≤ σ_SI ≤ max_x c(x).   (3.31)

**Angular vs spatial discretization, general verdict** (Ch. III summary, p. 69, items 3–4):

> "The discrete-ordinates approximation of the angular variable … generally has only a
> minor effect on convergence rates, and this effect does not lead to divergence for
> standard symmetric angular quadrature sets." / "For the first-order form … the choice
> of the spatial discretization schemes for the high-order and low-order DSA equations
> has a much more significant effect on convergence rates than the choice of the angular
> quadrature set."

### 2.3 DSA, continuous theory: the 0.2247c result

**Scheme** (§II.B, pp. 27–28): sweep (2.22)–(2.24) → exact correction equation
(2.27)–(2.29) → replace by P₁/diffusion (2.30)–(2.32) → additive update
φ^(ℓ+1) = φ^(ℓ+1/2) + F^(ℓ+1) (2.34). "Diffusion Preconditioning" is the review's
alternative name (p. 27).

**Fourier system** (Eqs. (2.43)–(2.46), p. 29, transcribed from the page image):

  (iλμ + 1)α(λ,μ) = c/2,
  β(λ) = ∫₋₁¹ α(λ,μ)dμ,
  (λ²/3 + 1 − c)γ(λ) = c[β(λ) − 1],
  ω(λ) = β(λ) + γ(λ).

**Closed forms** (Eqs. (2.47)–(2.50), p. 29):

  α(λ,μ) = (c/2)·1/(1+iλμ),
  β(λ) = (c/λ)tan⁻¹λ    [the SI eigenvalue],
  γ(λ) = 3c/(λ²+3(1−c)) · ((c/λ)tan⁻¹λ − 1),

  **ω(λ) = c·[λ²/(λ²+3(1−c))]·∫₀¹ (1−3μ²)/(1+λ²μ²) dμ
        = 3c/(λ²+3(1−c)) · [(λ²/3 + 1)·(tan⁻¹λ)/λ − 1].   (2.50)**

**The headline numbers** (p. 29, verbatim):

> "For c = 1, the spectral radius of the SI scheme σ_SI = 1, while the spectral radius of
> the DSA scheme σ_DSA ≈ 0.2247. For all λ, 0 ≤ ω(λ) ≤ β(λ), and ω → 0 as λ → 0 or
> λ → ∞. Thus, the DSA scheme efficiently suppresses the error modes that have strong
> spatial and angular dependence *and* those that have weak spatial and angular
> dependence."

  ω_c(λ) ≤ ω₁(λ)·c ≤ 0.2247c.   (2.51, p. 29)

Conditions for 0.2247c: infinite homogeneous medium, isotropic scattering, one group,
continuous (undiscretized) equations, P₁ low-order closure. The bound is attained at an
**intermediate** λ (Figure 1, p. 30, plots β and ω for c = 1: SI peaks at λ = 0, DSA peaks
near λ ≈ 2.5); the c-scaling (2.51) is "sharp only for c ≈ 1" (p. 30).

**Iteration-count corollaries** (pp. 30–31): 10^−m ≈ (0.2247c)^{ℓ_DSA} (2.52) ⇒

  ℓ_DSA ≤ [ln10/(ln(1/c) + ln(1/0.2247))]·m ≤ 1.54·m   (2.53)
  ℓ_SI = [ln10/ln(1/c)]·m   (2.54);  c = 0.999 ⇒ ℓ_SI = 2301·m.

**Cost model** (p. 31): in TWODANT a multi-D diffusion solve ≈ one S₄ sweep ⇒ one DSA
iteration ≈ 2 SI iterations ⇒ DSA beats SI for roughly 0.1 < c < 1 (with (2.50): "DSA is
much more efficient than SI for c ≈ 1 and remains more efficient for c as small as 0.1").

### 2.4 The consistency requirement — the load-bearing DSA story

This is the review's central technical narrative (§III.B, pp. 52–58 + history §I.G,
pp. 16–17). Sequence:

**(a) Inconsistent DSA behavior.** Discretize the transport with DD-S_N but the diffusion
with a conventional cell-edge scheme (Eqs. (3.32)–(3.34), p. 52). Fourier result
(Eqs. (3.39)–(3.41), p. 53, from the page image):

  ω = ω_SI − c(1−ω_SI)·cos²(λΣ_t h/2) / [1 − c + ⅓((2/(Σ_t h))sin(λΣ_t h/2))²]   (3.40)

  ω = ω_SI − c(1−ω_SI) / {(1−c)[1 + (Σ_t h/2)²Λ²] + ⅓Λ²}.   (3.41)

At fixed Λ, ω is **monotonically increasing in Σ_t h** (p. 53):

  ω_min(Λ) = ω_SI − c(1−ω_SI)/(1 − c + ⅓Λ²)  at Σ_t h → 0,   (3.43)
  ω_max(Λ) = ω_SI  at Σ_t h → ∞.   (3.44)

> "As Σ_t h → ∞, Eq. (3.44) holds and the performance of the DSA scheme reduces to
> that of SI. In this case, the discrete DSA scheme is convergent but ineffective — it
> performs like Source Iteration." (p. 53; Figure 5, p. 54: σ_DSA vs Σ_t h for S₂/S₄/S₈,
> c = 0.98.)

**(b) Worse: instability for other common discretizations** (p. 54, verbatim — the M2
anchor):

> "This result is marginally acceptable, but much worse results hold for many other
> common discretizations of the low-order DSA diffusion problem: the resulting DSA
> scheme is rapidly convergent for small Σ_t h, but it degrades as Σ_t h increases and then
> becomes divergent or unstable (σ_DSA > 1) for Σ_t h sufficiently large. This instability
> occurs, for example, with common cell-centered diffusion discretizations [51]."

([51] = Reed 1971, NSE 45:245 — the original demonstration; history p. 16: "Reed then
showed by testing and analysis that while the 'diffusion' synthetic method employed by
Gelbard and Hageman is rapidly convergent for sufficiently fine spatial grids, it diverges
for coarser grids.")

**(c) Alcouffe's remedy and the Four-Step derivation** (pp. 54–57). Alcouffe [57, 63]
solved it for DD; the review presents the general recipe as the (Larsen) "Four-Step"
procedure, Eqs. (3.45)–(3.64):

1. Formulate the **exact discrete** S_N problem for the corrections f (Eqs. (3.46)–(3.49)).
2. Operate on it with the **discrete P₁ approximation**: take the w_n and μ_n w_n moments
   and insert f_{n,j(±1/2)} ≈ ½(F_{j(±1/2)} + 3μ_n G_{j(±1/2)}) (3.50). Boundary detail:
   the discrete Marshak-like condition carries γ ≡ Σ_{μ_n>0} μ_n w_n → ½ only as N → ∞
   (3.54), (3.56).
3. Algebraically eliminate everything except the cell-edge scalar-flux corrections →
   a tridiagonal edge-centered diffusion system (3.61)–(3.63). Its differences from the
   "conventional" discretization (p. 56): a **three-point** (not one-point) removal (Σ_a)
   term, and the γ_N boundary factor.
4. Solve and update (3.64).

**(d) The consistent result** (p. 57):

  ω = ω_SI − c(1−ω_SI)/(1 − c + ⅓Λ²)  **for all Σ_t h**   (3.65)

> "This result is identical to Eq. (3.43), which we have shown represents a
> rapidly-convergent scheme! … The 'consistent' spectral radii are independent of
> Σ_t h, and for all quadrature orders N and all 0 < c ≤ 1.0 satisfy σ_DSA < 0.2247."
> (p. 57; Figure 6 compares inconsistent vs consistent σ_DSA vs Σ_t h.)

Interpretation (p. 57): "a way to obtain a low-order DSA diffusion discretization that
accelerates effectively for any Σ_t h is to **derive this discretization from the discrete
S_N equations, using the discrete P₁ approximation**." The Four-Step method is also
identified as "a direct application of Lebedev's KP1(1) method to the fully-discrete
transport problem" (p. 58), with the operator decomposition M = E D⁻¹K of Eqs.
(1.37)–(1.40) (K = P₁ projection, D = the low-order operator obtained *by applying K to
L−S*, E = prolongation).

**(e) The canonical WHY-consistency paragraph** (p. 58, verbatim — quote this in the docs):

> "In DSA, the importance of a consistently-discretized low-order diffusion equation
> cannot be overstated. Generally, DSA with a conventionally-discretized low-order
> diffusion equation becomes inefficient or even divergent as Σ_t h increases. To
> maintain rapid convergence for all Σ_t h, careful attention must be paid to the
> low-order diffusion discretization."

Practical performance of consistent DSA (p. 58): "the iterative error reduces at least two
orders of magnitude for every three DSA iterations. … most practical 1-D problems can be
solved with a computational cost less than that of about 12 source iterations."

**(f) The operator-level cause** (Ch. I, p. 21 — why the low-order operator is a knife
edge): in the preconditioned iteration matrix X_SI − MS(I−X_SI),

> "I − X_SI has small eigenvalues when X_SI has eigenvalues close to unity. Thus, for
> such error modes, MS must have large eigenvalues to keep the iteration matrix … from
> approaching the unaccelerated matrix X_SI. However, if the eigenvalues of MS are too
> large, then the scheme diverges. Thus, the operator M must be chosen very carefully."

The review names the resulting engineering burden "the **tyranny of consistent low-order
discretization**" (p. 20).

**(g) Escapes from the tyranny** (each with its price):

- **S₂SA** (§II.C pp. 31–32, §III.C p. 58): low-order = S₂ with the *same* spatial
  scheme ⇒ "eliminates the issue of consistent diffusion differencing" (p. 32); in 1-D
  algebraically ≡ DSA (P₁ ≡ S₂ in slab). Price: costlier low-order solve; **not**
  equivalent in 2-D/3-D (σ_S2SA → c/(2−c), see 2.6).
- **TSA** (§II.E pp. 34–37): low-order = transport with Σ_s → (1−β)Σ_s, same
  discretization ⇒ "essentially independent of the spatial discretization scheme" (p. 77).
  Price: inner iterations; β trade-off (see 2.6).
- **Nonlinear/QD-family** (see 2.8): rapid convergence *without* consistency, but the
  converged answer is no longer the unaccelerated S_N solution.
- **Second-order (even-parity/SAAF) forms** (§II.I pp. 45–48, §III.I, §IV.E): "Because
  the second-order transport equations are already in diffusion-like form, it is
  straightforward to discretize the diffusion equation consistently with them" (p. 48);
  Lebedev's KP success lived here (p. 141). Price: no sweep solve ("cannot be solved by
  simple transport sweeps," p. 48), banded systems, void-region trouble (D = 1/(3Σ_t) → ∞).

### 2.5 Multi-D theory: rotation argument, quadrature-orientation effects, and where 0.2247c survives

**Continuous, any dimension** (§IV.D.1, pp. 80–81): each Fourier error mode is a plane
wave; a coordinate rotation reduces it to slab ⇒ SI and DSA behave **identically to 1-D**
in continuous multi-D Cartesian geometry; "The behavior of an error mode depends only
upon its wave number λ … not on its orientation" (p. 81).

**Discrete ordinates break rotational invariance** (p. 81):

> "now one cannot perform a coordinate rotation … The reason for this is that the
> discrete-ordinates operator is not rotationally invariant. Thus, SI and DSA may
> perform differently for error modes that have the same wave number but different
> orientations. … we would expect the largest differences to occur with quadrature sets
> that have few directions."

- **SI multi-D** (Eq. (4.59), p. 82): ω(λ⃗) = c Σ_n w_n/(1+(λ⃗·Ω̂_n)²)·(normalized); flat
  mode still gives c. Orientation effect quantified for S₂ (Eqs. (4.60)–(4.61)):
  mode along a quadrature direction ~7c/9λ² vs along a coordinate axis ~3c/9λ² at large λ.
- **DSA multi-D** (Eq. (4.65), p. 83): ω = [3c/(λ²+3(1−c))]·Σ_n w_n[λ² − 3(λ⃗·Ω̂_n)²]/[1+(λ⃗·Ω̂_n)²]·(norm).
  For S₂: ω = 0 for axis-aligned modes (4.66) but **ω → c/2** as λ → ∞ for diagonal
  modes (4.67). Verbatim (p. 83):

  > "We recall that in 1-D the DSA spectral radius is zero for the S₂ quadrature set and
  > monotonically increases with increasing quadrature resolution up to approximately
  > 0.2247c for S∞. The situation in 2-D and 3-D is entirely different: the DSA spectral
  > radius is 0.5c for the S₂ quadrature set and monotonically decreases with increasing
  > quadrature order down to approximately 0.2247c for S∞." (Detailed study: Adams &
  > Wareing [181].)

- **S₂SA multi-D** (Eq. (4.68), p. 83, Gelbard-Hageman [32]): σ_S2SA → c/(2−c) as the
  high-order quadrature refines (→ 1 as c → 1); S₄SA: c/(6−c) < 0.2c (Larsen [94]);
  two sweeps per S₂ solve (Lebedev K₂P): ≤ 0.153c (p. 84).
- **Consistent multi-D DSA exists for DD**: the four-step vertex-centered 9-point (2-D) /
  27-point (3-D) diffusion discretization "yields a DSA spectral radius equal to that
  without spatial discretization" (p. 85); production codes (DANTSYS/TWODANT) instead use
  a **slightly inconsistent 5/7-point** approximation to cut cost — "a winning strategy
  for problems that are not too optically thick and highly scattering, but it does not
  work efficiently for the most difficult problems" (p. 87).
- For non-DD multi-D discretizations (BLD/corner-balance): full four-step is
  intractable-to-solve (12 unknowns/cell, cannot eliminate, p. 85); the modified-DSA
  family (Khalil; Adams-Martin M4S; Wareing symmetric-simplified; asymptotic DSA) trades
  a little consistency for solvability — see Task 3 for their failure modes. First
  **unconditionally efficient** non-DD multi-D scheme: Morel-Dendy-Wareing 1993 [190]
  (BLD transport → BLD diffusion line relaxations → BLC diffusion multigrid; "converges
  any fixed-source one-group isotropic-scattering transport problem for the cost of a few
  transport sweeps," p. 88).
- **Unstructured grids** (§IV.D.5, pp. 88–89): open. "at the time of this writing, an
  unconditionally efficient acceleration scheme for transport on unstructured grids has
  not been demonstrated" (p. 89).

### 2.6 Taxonomy of iteration/acceleration families (as the review organizes it)

Linear preconditioning ("synthetic") family — sweep + additive low-order correction;
converged answer = unaccelerated S_N solution; consistency required:

| Method | Low-order operator | Key spectral result | Where |
|---|---|---|---|
| SI (Richardson) | none | ρ = c (any N, any h) | §II.A, §III.A |
| DSA = KP1(1) | P₁/diffusion | ≤ 0.2247c (consistent); Σ_t h-degradation or divergence (inconsistent) | §II.B, §III.B |
| S₂SA | S₂, same spatial scheme | ≡ DSA in 1-D; c/(2−c) in multi-D | §II.C, §III.C, §IV.D.3 |
| KP1(1) w/ tuned D-coefficient | q/Σ_t, q free | q = 1/3: ≤ 0.225c; q = 0.281: ≤ 0.186c (Lebedev; Azmy [160] similar) | §II.D p. 32 |
| KP1(1)P₂(0) | diffusion + algebraic P₂ stage | ≤ 0.127c (Lebedev) | §I.F p. 15, §II.D |
| TSA | transport w/ Σ_s → (1−β)Σ_s | inner ρ = (1−β)c/(1−βc) (2.67); β = 1 ⇒ Two-Step (Morozov [28]/Larsen-Miller [110]); outer ρ grows with β; M-sweep truncation (2.73) | §II.E, §III.E, §IV.A.3 |
| BPA / ICSA | angular projection **on cell surfaces only** (double-P₀ etc.) | projection-operator-dependent: c/3, 0.2247c, or divergent (p. 96) | §V.B |
| Upscatter Two-Grid | 1-group diffusion w/ spectral-shape weighting ξ_g | heavy water: GS 0.9998 → TG 0.489 (= 2nd GS eigenvalue) | §VII.A |
| Angular multigrid (forward-peaked) | recursive S_{N/2} + DSA base | ≤ 0.6 (1-D, Fokker-Planck limit); multi-D needs diffusion filtering of divergent high-frequency modes | §VII.B |

Nonlinear family — low-order equation solves directly for the accelerated quantity
(multiplicative/functional closure); no consistency needed; answer ≠ unaccelerated S_N
solution unless consistent; natural for k-eigenvalue:

| Method | Closure functional | Notes | Where |
|---|---|---|---|
| Quasidiffusion (QD, Gol'din 1964 [25]) | Eddington factor E = ∫μ²ψ/∫ψ | linearized-QD ≈ Lewis-Miller/DSA (p. 40); low-order eigenvalue-ready; nonsymmetric low-order | §II.G, §III.G, §VIII.D |
| Lewis-Miller (LM) [61] | linear, non-preconditioning (low-order directly yields next iterate) | linear-but-not-true-acceleration; Eq. (1.46) shows any D-discretization converges | §II.F, §III.F |
| Weighted-Alpha (WA) family (Anistratov-Larsen) | \|μ\|^α half-range functionals | σ: α=0 (First-Flux) 0.3333; α=1 (Second-Flux) 0.3105; α≈0.366 (S₂-like) 0.2247; α≈0.128 optimal 0.1865 (p. 45) | §II.H, §III.H |
| Nonlinear BPA (Adams [178]) | boundary angular-shape functionals | linearizes to BPA; high/low-order discretizations may differ deliberately (accuracy-enhancement use, Fig. 12 p. 97) | §V.B |
| Nonlinear rebalance / CMR | multiplicative per-coarse-cell factors | see Task 3 — divergence windows | §V.A |
| Alcouffe's nonlinear k-DSA | D and Σ_a rational functionals of transport solution | DANTSYS production scheme; fixup interactions (Task 3) | §VIII.E |

Algebraic family — no low-order physics operator; Ch. VI:

- **ASE/Lyusternik** (§VI.A p. 103): single-mode extrapolation; "can be divergent … if
  ASE is applied after every source iteration in a problem with c > 1/2"; "not used in
  modern S_N codes."
- **Chebyshev** (§VI.B pp. 104–107): polynomial-averaged Richardson iterates;
  σ_Cheby = (√κ−1)/(√κ+1) (6.19); needs eigenvalue-interval estimates (a, b);
  two-term cyclic variant "suffers from severe numerical instabilities" unless the
  linear factors are carefully ordered (Lebedev-Finogenov / Samarsky-Nikolaev, p. 107).
- **CG** (§VI.C pp. 107–111): same asymptotic rate as optimally-tuned Chebyshev without
  eigenvalue estimates; requires SPD ⇒ isotropic scattering + symmetric quadrature; the
  transport-correct inner product carries a Σ_s weight, ⟨f,g⟩ = ∫_{Σs>0} Σ_t c f g d³r
  (6.32) (Sanchez-Santandrea symmetrization [298] noted as the anisotropic escape).
- **Krylov (GMRES/BCGSTAB)** (§VI.D pp. 111–113): see 2.7.
- The unifying identity (§VI.B p. 106): **κ ≈ 1/(1−σ_SI)** (6.20) —

  > "the statements 'SI is slow to converge' and 'T has a large condition number' are
  > equivalent. … The way to make Chebyshev acceleration fast is to find a good
  > preconditioner … this is the same problem of finding a good low-order operator in a
  > 'synthetic' acceleration scheme — the problem that has been the subject of much of
  > this paper." (pp. 106)

Other iterative structures: rebalance (§V.A), spatial multigrid (§V.C, pp. 98–103:
excellent in 1-D for first- and second-order forms; 2-D extensions in progress),
k-eigenvalue schemes (Ch. VIII; see Task 1 row).

### 2.7 Krylov + DSA-as-preconditioner — the ORPHEUS-relevant results

**The linear-algebra identification** (pp. 12–13, Eqs. (1.27)–(1.34)): SI ≡ Richardson
iteration on Aψ = η̂, A ≡ I − L⁻¹S, η̂ = L⁻¹q; synthetic acceleration ≡ **preconditioned
Richardson** with P = I + MS ≈ A⁻¹ (1.32)–(1.34). Scalar-flux form (p. 21, Eq. (1.45)).
Chapter VI operator (p. 104, Eqs. (6.4)–(6.5)): T = I − X_SI with X_SI = K₀L⁻¹(Σ_s/2)·(·);
Tφ = η̂, η̂ = K₀L⁻¹q — **this is exactly the ORPHEUS GMRES-on-the-honest-matvec system
after sweep preconditioning** (one sweep = one application of X_SI; T never formed:
"we simply execute a transport sweep … and then perform a simple manipulation," p. 107;
the matrix-free matvec is Td = d − K₀L⁻¹(Σ_s/2)d, Eq. (6.26)).

History (pp. 110–112): Derstine & Gelbard [100] first (1985, CG+DSA); Faber & Manteuffel
[139] independently recognized DSA = preconditioned Richardson; Ashby, Brown, Dorr &
Hindmarsh [149, 158, 184, 204] did the systematic theory+experiments; Ashby et al. were
"the first to apply Krylov methods to transport problems with anisotropic scattering
(which causes T to be asymmetric)" via GMRES (p. 112).

**Finding 1 — Krylov rescues inconsistent DSA** (Derstine-Gelbard Table 1, p. 110;
XY geometry, DD transport, cell-centered FD diffusion, c = 0.99 in part of the problem):

| Σ_tΔx | 0.5 | 1.0 | 2.0 | 4.0 |
|---|---|---|---|---|
| TWODANT DSA | 6 | 6 | 7 | 12 |
| Inconsistent DSA + CG | 10 | 12 | 13 | 18 |
| Inconsistent DSA alone | 12 | 15 | diverged | diverged |

> "for some problems, inconsistent DSA diverges by itself but with CG it converges in
> less than 20 iterations!" (p. 110). Mechanism: CG/GMRES convergence is "theoretically
> guaranteed … given SPD matrices" regardless of preconditioner quality — wrapping the
> unstable fixed-point map in a Krylov solver converts divergent modes into merely
> slow ones.

**Finding 2 — Krylov does not replace a good preconditioner** (Ashby et al., pp. 110–113).
Verbatim (p. 110–111):

> "perhaps most importantly, a 'good' DSA method is much faster to converge by itself
> than is a 'poor' method with CG. Further, because a 'good' DSA scheme (such as one
> arising from Larsen's Four-Step Method) converges so quickly, employing CG in
> addition generally reduces iteration counts by only one or two. Thus, it seems that
> while CG does help, it does not remove the need for good preconditioners."

Table 2 (p. 113, "Ashby's problem": two-region slab, 50 mfp c=0.99998 + 5×10⁷ mfp c=0.8,
16-pt Gauss, tol 10⁻⁸ on residual, cap 50):

| Transport | Diffusion precond. | Richardson | Chebyshev | GMRES |
|---|---|---|---|---|
| DD | consistent | 11 | 8 | 7 |
| DD | FEM (inconsistent) | >50 | 22 | 13 |
| LD | consistent | 11 | 8 | 7 |
| LD | FEM (inconsistent) | 43 | 20 | 14 |

Second Ashby problem (uniform slab, 10⁸ mfp, c = 0.998): with the FEM preconditioner
**all** methods failed at 50 iterations; with the consistent preconditioner all converged
in ≤ 2 (p. 113). Verdict verbatim (p. 113):

> "while GMRES and Chebyshev are improvements over Richardson iteration, nothing takes
> the place of a good preconditioner (such as a consistent diffusion discretization)."

**Finding 3 — the design guidance** (Discussion, p. 140, verbatim — this is the sentence
that licenses the ORPHEUS architecture of GMRES-with-sweep now (#200) + DSA later (#2)):

> "A trend that has received increasing attention for difficult problems — ones with,
> perhaps, severe spatial discontinuities, unstructured spatial grids, or even a slightly
> inconsistent acceleration scheme — is to employ an algebraic iterative method such as
> Conjugate Gradient (CG) or GMRES 'outside of' DSA (i.e., to use CG or GMRES with a
> diffusion preconditioner). If the sweep + preconditioner has only 'a few'
> slowly-converging or unstable error modes, then this approach can be quite
> advantageous. However, it remains true that for maximum overall efficiency, the
> diffusion preconditioner should suppress the maximum number of slowly converging
> error modes, without causing other modes to become slowly-converging or unstable, at
> the minimum cost. (If the diffusion preconditioner efficiently suppresses all error
> modes, then it is unnecessary — but probably not harmful — to combine it with an
> algebraic method such as CG or GMRES.)"

Secondary Krylov use (p. 113): Krylov **inside** the preconditioning step (solving an
asymmetric low-order problem, e.g. Hong-Cho NLBPA low-order via BCGSTAB [252]); Ramoné's
TSA used CG to solve the low-order transport problem rather than to accelerate the outer
iteration (p. 111); Zika-Adams TSA with opposing reflecting BCs recovered CG-compatibility
by identifying the underlying symmetric infinite-medium problem [284, 285] (p. 111).

### 2.8 Nonlinear methods (CMFD/NDA lineage)

The review predates the "CMFD"/"NDA" labels but contains the lineage as **Quasidiffusion
(Gol'din 1964 [25])**, the **flux/WA methods**, **nonlinear rebalance**, **NLBPA**, and
**Alcouffe's nonlinear k-DSA** (the D̂-style transport-corrected diffusion functional of
§VIII.E is the direct ancestor of CMFD-class closures; Smith-style CMFD itself is not
cited — the review's coarse-mesh nonlinear content is rebalance + QD).

The linear-vs-nonlinear contract (Discussion pp. 138–139, condensed from the two bullet
lists — docs-ready):

- Linear schemes: additive correction; converged answer = unaccelerated S_N solution;
  **"for unconditional convergence, require the low-order preconditioning equation to be
  discretized consistently"**; do not require positivity; "natural to apply to
  source-detector problems, but are not natural for eigenvalue problems."
- Nonlinear schemes: low-order equation for the solution itself; **two** converged
  estimates differing by truncation error; **no consistency requirement for rapid
  convergence** (but consistent discretization ⇒ both estimates = unaccelerated
  solution); "generally require the converged solution and all iterates leading to this
  solution to be positive"; "natural to apply to both source-detector and eigenvalue
  problems."

QD convergence detail (§III.G, p. 65 + summary item 7, p. 70): rapidly convergent for all
spatial grids without consistency, but "unlike preconditioning methods, they are not true
acceleration methods; upon iterative convergence, they do not yield the unaccelerated
solution of the spatially discrete S_N equations" — framed as **either contamination or
opportunity** ("a good low-order discretization can produce improvements over the
unaccelerated solution"; streaming/void problems → contamination, diffusive problems →
opportunity, p. 70). The QD philosophy, stated at the k-eigenvalue application (p. 132):
"most of the work is done in the relatively inexpensive setting of a lower-dimensional
problem." Historical scoop (p. 141): the Soviet groups (Gol'din QD; Lebedev even-parity
KP) "avoided the difficulties of 'consistent discretization'" by construction — the
east/west split of the field is precisely over where the consistency burden lands.

### 2.9 Geometry-generality ledger (which convergence results hold where)

| Result | Slab 1-D | Sphere 1-D | Cylinder 1-D | Multi-D Cartesian |
|---|---|---|---|---|
| ρ_SI = c (continuous) | (2.17) | **identical** — proved via integral-equation trick, Eq. (4.39) p. 78 | not analyzed (no result; expected by analogy) | identical (plane-wave rotation, Eq. (4.53) p. 81) |
| ρ_DSA ≤ 0.2247c (continuous) | (2.50)–(2.51) | **identical**, Eq. (4.43) p. 79 — "To our knowledge, this analytic result is new" | not analyzed | identical (rotation argument, p. 81) |
| ρ_SI = c (discrete) | (3.30), any N, h | **no Fourier analysis possible** with the α-recursion angular discretization (p. 79) | none | flat mode still c; orientation-dependent otherwise (4.59)–(4.61) |
| ρ_DSA (discrete) | 0.2247c (consistent) | none — "development … has proceeded to date without the benefit of a sharp Fourier analysis" (p. 79) | none | S₂: c/2 → S∞: 0.2247c (4.65)–(4.67); heterogeneity/aspect-ratio degradation (Task 3) |

The spherical trick (docs-worthy method): the separable Fourier ansatz fails ("the error
is not separable. Thus, conventional wisdom has been that it is not possible to
Fourier-analyze iteration schemes in spherical geometry," p. 78) — instead eliminate the
angle analytically to get the integral equation for the scalar-flux error (Eq. (4.35),
E₁ kernel), extend rφ(r) oddly to r < 0, and use the ansatz rφ^(ℓ)(r) = ω^ℓ sin(λΣ_t r)
(4.38) ⇒ ω(λ) = (c/λ)tan⁻¹λ (4.39), identical to slab. Same route for DSA ⇒ (4.43) ≡
(2.50). Both results are **continuous-geometry only**; discretized-in-angle curvilinear
S_N (Eqs. (4.44)–(4.46), the standard α_{n±1/2} recursion) admits no known Fourier
analysis (p. 79). Curvilinear practice is empirical (p. 80): "Experience has been that the
techniques that work in slab geometry also work in spherical geometry. Thus, while
spherical geometry makes it difficult to analyze iterative methods, it does not seem to
hamper their effectiveness [144, 175]."

---

## Task 3 — M2 harvest: symptom → cause → fix triples

1. **False convergence** (§I.C, pp. 9–10) — *the universal one; belongs at the top of any
   ORPHEUS convergence-diagnostics section.*
   - Symptom: successive-iterate test ‖ψ^(ℓ)−ψ^(ℓ−1)‖ < ε passes, but the answer is wrong.
   - Cause: Eq. (1.19): ‖ψ−ψ^(ℓ)‖ ≤ ε·σ/(1−σ); for σ = 0.999 "the iteration error may be
     three orders of magnitude greater than the convergence criterion" (p. 10; ~2303
     iterations per decade, Eq. (1.20)).
   - Fix (p. 10): test against ε(1−σ̂) with σ̂ estimated on the fly,
     σ̂ ≈ ‖ψ^(ℓ+1)−ψ^(ℓ)‖/‖ψ^(ℓ)−ψ^(ℓ−1)‖ (valid when the dominant eigenvalue is real and
     simple, "usually true in transport problems").

2. **Inconsistent-DSA mesh-thickness degradation** (§III.B, p. 53).
   - Symptom: DSA fast on fine meshes, stalls toward SI-like counts as cells thicken.
   - Cause: Eq. (3.41) — ω monotonically increasing in Σ_t h; ω_max = ω_SI (3.44):
     "convergent but ineffective — it performs like Source Iteration."
   - Fix: consistent (four-step-derived) low-order discretization ⇒ (3.65),
     Σ_t h-independent, σ < 0.2247. User-facing alternative (Discussion, p. 138): if the
     grid is known never to be optically thick, a cheap inconsistent DSA is a legitimate
     design point — "the code user must be aware that if the code is applied to problems
     with increasingly thick spatial cells, the rate of iterative convergence of the code
     will degrade" (p. 71).

3. **Inconsistent-DSA divergence** (§III.B, p. 54; history p. 16). *This is the diagnosis
   Stacey's one-liner ("the synthetic method may even become unstable with small mesh
   spacing") lacks — note the review locates instability at COARSE mesh (Σ_t h large),
   not small; the Gelbard-Hageman/Reed story is precisely "rapidly convergent for
   sufficiently fine spatial grids, … diverges for coarser grids" (p. 16).*
   - Symptom: iteration diverges (σ_DSA > 1) when cells are optically thick, worse as c→1.
   - Cause: the low-order diffusion operator over-corrects modes the discrete sweep
     attenuates differently than the analytic sweep — "common cell-centered diffusion
     discretizations [51]" (Reed 1971).
   - Fix: (i) consistency (Alcouffe [57, 63]; four-step [83, 85]); (ii) same-discretization
     transport low-order (S₂SA/TSA, §III.C p. 58); (iii) wrap in Krylov — rescues
     convergence but not efficiency (Table 1 p. 110, Table 2 p. 113); (iv) nonlinear
     closure (QD/WA), which "do not generally experience iterative divergence" (p. 70).

4. **Consistent-but-lagged variants: γ_N boundary factor** (p. 57). Minor but real: the
   consistent scheme's boundary condition uses γ_N = Σ_{μ_n>0}μ_n w_n, which "limits to
   the Eq. (3.33) value of 1/2 only as N → ∞" — using ½ (the continuum Marshak value) in
   place of γ_N is a quiet consistency break at small N.

5. **Multi-D DSA degradation from heterogeneity** (§III.J p. 69; Discussion p. 139; refs
   [259], [260], [279] — the Azmy adjacent-cell line).
   - Symptom: even *fully consistent* DSA converges slowly (observed σ up to 0.88) on 2-D/3-D
     problems.
   - Cause: "characterized as being 2-D or 3-D (not 1-D), with strong spatial
     heterogeneities, and with spatial grids having large aspect ratios. It is believed
     that for such problems, the diffusion approximation to the transport equation is not
     sufficiently accurate to provide an efficient preconditioner" (p. 69). Discussion
     p. 139: "spectral radii can grow as large as 0.88, even with methods that, for
     homogeneous problems, always produce spectral radii smaller than 1/3. **We are not
     aware of a remedy to this interesting degradation.**"
   - Status marker (p. 69): "(We are not aware that any problems have been seen in which
     a 'consistently-discretized' DSA scheme has actually become unstable.)" — i.e.,
     degradation yes, instability no (as of 2002). [Post-review literature (Warsa-Wareing-
     Morel 2004) continues this thread; outside this survey's source.]

6. **Asymptotic-DSA aspect-ratio failure** (§IV.D.4, pp. 86–87).
   - Symptom: Wareing-Larsen-Adams "asymptotic DSA" effective in 1-D and on 2-D square
     cells, stalls on stretched cells.
   - Cause: "as the 2-D cell aspect ratio approaches infinity, the spectral radius of
     asymptotic DSA approaches c" (p. 87); root mechanism (p. 86): "error modes that are
     discontinuous at cell edges; the diffusion-limit spatial discretization does not
     attenuate these modes, and their eigenvalues approach c for optically thick cells."
   - Fix: local (cell-independent) post-diffusion discontinuous correction (p. 86–87) —
     fixes the thick-cell limit but not the aspect-ratio limit.

7. **M4S divergence on skewed grids** (§IV.D.5, p. 88).
   - Symptom: Adams-Martin M4S (rapidly convergent on orthogonal grids) diverges.
   - Cause: skewed tetrahedral cells (Warsa et al. [300]).
   - Moral (p. 88): "one should not assume that the orthogonal-grid behavior of an
     iterative scheme will automatically extend to arbitrary grids." Companion cost
     finding (p. 89): fully consistent four-step on tets has small σ but is "extremely
     expensive relative to a transport sweep" — small spectral radius and affordable
     solve are separate requirements (restated p. 139).

8. **Rebalance divergence windows** (§V.A, pp. 90–95; Figures 9–11; linearized Fourier
   analysis of Cefus-Larsen [127, 150]). Five numbered facts (p. 93):
   - "Rebalance is generally divergent if the optical thickness of the spatial cells Σ_t h
     is too small. In particular, rebalance is very difficult to employ effectively if the
     problem contains voided regions."
   - "usually convergent for scattering ratios c sufficiently small. However, as c → 1,
     rebalance often becomes divergent. Unfortunately, these are the problems for which
     effective acceleration is most needed."
   - "Fine-mesh rebalance (p = 1) is usually divergent for c ≈ 1, unless Σ_t h ≫ 1."
   - Coarse-mesh: convergent for enough fine cells per coarse cell (large p); too large p
     ⇒ "convergent but inefficient (the convergence rate converges to that of source
     iteration)"; optimal p is a user-tuned window between divergence and uselessness.
   - Diagnosed cause: "The coarse-mesh rebalance equations do not limit to a physically
     or mathematically meaningful result as Σ_t h → 0. (This is likely why rebalance
     becomes divergent in this limit.)"
   - Fix: none within rebalance — "the particle transport community has generally come to
     prefer iterative acceleration methods such as DSA" (p. 95).

9. **TSA divergence at large β** (§II.E, p. 35).
   - Symptom: outer iteration diverges as the low-order scattering reduction β → 1.
   - Cause: "For more complicated problems in multidimensional geometries or with
     anisotropic scattering, the spectral radius can exceed unity for β close to unity"
     [243]; and with fixed inner-sweep count M, "for optically thick problems with very
     small absorption, the spectral radius of the TSA scheme (with fixed M) still
     approaches unity" (p. 35).
   - Fix: keep β small enough (problem-dependent trade against inner cost, Eq. (2.67)).

10. **Negative-flux fixups × acceleration** (§VIII.E, pp. 136–137) — the review's explicit
    fixup-interaction diagnostic.
    - Symptom: accelerated eigenvalue iteration degrades or code fails when fixups engage.
    - Cause: nonlinear low-order functionals "require positive transport solutions, which
      almost always translates into a requirement for a negative-flux fixup … First, in
      eigenvalue problems, negative flux fixups tend to decrease solution accuracy.
      Second, **fixups can cause a transport discretization to be inconsistent with a
      diffusion discretization that was carefully constructed to be consistent**" (p. 137).
      Related in-family failure: Alcouffe's D-functional (8.46) "can produce a negative or
      very large diffusion coefficient" ⇒ fallback removal-correction (8.47); DANTSYS
      "permits" negative corrected Σ_a "unless it causes the discrete diffusion matrix to
      lose diagonal dominance. In such cases … the code simply turns off the acceleration
      in the offending energy group" [287] (p. 136).
    - Fix: linear k-DSA (Gelbard et al. [81]) — "remained robust and fast (without
      negative flux fixups), even on problems for which negative fixups caused the DANT DD
      code to fail"; costs one adjoint low-order eigenvalue problem + a singular-operator
      low-order solve per iteration (p. 137).

11. **Gauss-Seidel upscatter stagnation** (§VII.A, pp. 114–117).
    - Symptom: multigroup thermal iterations stall in large weakly-absorbing thermal systems.
    - Cause: with zero absorption the GS iteration has an eigenvalue → 1 at λ = 0 whose
      eigenvector "is a discrete approximation of a Maxwellian" (p. 115) — the physical
      thermalization equilibrium is the slow mode (heavy water 41-group: ρ_GS = 0.9998).
    - Fix: Adams-Morel Two-Grid [179, 180] (one-group diffusion correction weighted by the
      slowest mode's spectrum ξ_g, precomputed per material from the λ = 0 eigenproblem);
      result ρ_TG = 0.489 for heavy water = the *second* GS eigenvalue (p. 117);
      Averin-Voloschenko variant [193]. Consistency carries over: "if one knows how to
      discretize DSA for a within-group scattering problem for a given transport
      discretization, one can apply the same discretization" to the TG equation (p. 117).

12. **DSA ineffectiveness under forward-peaked scattering** (§VII.B, pp. 118–121).
    - Symptom: DSA barely better than SI for electron/optical-photon transport.
    - Cause: mode-by-mode SI eigenvalues at λ = 0 are ρ_n = Σ_sn/Σ_t (Eq. (7.14)); the
      diffusion preconditioner "addresses only the φ₀ and φ₁ portions of the scattering
      source. Thus, if c_n is close to unity for n > 1, even DSA will not cause the
      iteration to converge quickly" (p. 119). Sharper multi-D hazard: "P₁ acceleration of
      S₄ calculations can be divergent in multidimensional problems for sufficiently
      highly forward-peaked scattering" (Adams-Wareing [181]; Marchuk-Lebedev [10])
      (p. 120); Morel's within-group Fokker-Planck-type coupling is the 1982 [86] context
      (see Task 5).
    - Fix: Manteuffel-Morel angular multigrid [164] (σ ≲ 0.6 worst-case 1-D); in multi-D
      the naive extension has "divergent high-frequency modes … the multigrid method
      amplifies these error modes" (p. 120–121) — Pautz et al. [270] fix by diffusion
      *filtering* of the corrections (Eq. (7.16)) plus Galerkin quadrature [143, 269].

13. **BPA projection sensitivity** (§V.B, p. 96) — a warning ORPHEUS's interface-current
    experience independently corroborates: with the same double-P₀ trial space,
    "different projection operators produced spectral radii of c/3, 0.2247c, and a value
    larger than unity (divergence)." The projection (test) choice, not just the subspace,
    decides convergence.

14. **Chebyshev parameter fragility** (§VI.B, pp. 106–107). Performance "depends
    sensitively on the value chosen for a; it needs to be as large as possible while not
    violating a ≤ λ_min"; and the two-term cyclic implementation is roundoff-unstable
    without the Lebedev-Finogenov/Samarsky-Nikolaev factor ordering (p. 107). ASE:
    divergent if applied every iteration with c > 1/2 (p. 103).

---

## Task 4 — Section structure and scope

**Scope statement** (p. 6): deterministic S_N problems only; the review is about step 2 of
the two-step deterministic paradigm ("The resulting discrete system of equations is
solved… iterative methods must be used. The best iterative methods are unconditionally and
rapidly convergent."); excludes P_N/spherical-harmonics iteration ("different mathematical
properties"). Aim (p. 23): "to describe and illustrate the main concepts that have led to
important practical results," not encyclopedic coverage; Chapters II–III are explicitly
"tutorial in nature." Completeness claim (p. 140): "To our knowledge, the list of such
[English-language] references given below is complete, up to the time of the writing of
this Review," plus a deliberate program of integrating the Russian literature (Lebedev,
Gol'din; assisted by Anistratov — Acknowledgements p. 141).

**Full outline** (Table of Contents, pp. 4–5; page ranges as printed):

- **I. INTRODUCTION** (6–23): A Overview 6 · B Source Iteration 8 · C False Convergence 9
  · D Early Acceleration (Chebyshev, Rebalance) 11 · E Preconditioning ("KP"/"Synthetic")
  Methods 11 · F "KP" Methods — Historical Review 14 · G "Synthetic" Methods — Historical
  Review 16 · H Quasidiffusion and Related — Historical Review 19 · I Multigrid Methods 22
  · J Review Outline 23
- **II. ITERATION SCHEMES FOR CONTINUOUS TRANSPORT PROBLEMS IN PLANAR GEOMETRY** (23–48):
  A SI 23 · B DSA 27 · C S₂SA 31 · D KP Methods 32 · E TSA 34 · F Lewis-Miller 37 ·
  G Quasidiffusion 40 · H Weighted-Alpha 43 · I Second-Order Forms 45
- **III. THE EFFECT OF SPATIAL AND ANGULAR DISCRETIZATIONS IN PLANAR GEOMETRY** (48–71):
  A SI 48 · B DSA 52 · C S₂SA 58 · D KP 59 · E TSA 62 · F LM 63 · G QD 65 · H WA 66 ·
  I Second-Order Forms 66 · J Performance in Non-Model Problems 68 · K Summary 69
- **IV. OTHER DIFFERENCING SCHEMES AND GEOMETRIES** (71–90): A 1-D Planar (Four-Step 73;
  Other DSA 74; TSA 77; LM/QD/WA 77) · B 1-D Spherical 77 · C 1-D Cylindrical 80 ·
  D Multidimensional Cartesian (SI 80; DSA 82; S₂SA 83; Spatial Discretization: Regular
  Grids 85; Unstructured Grids 88) · E Second-Order Forms in Multi-D 89
- **V. OTHER ITERATIVE METHODS** (90–103): A Rebalance 90 · B Boundary Projection
  Acceleration 95 · C Spatial Multigrid 98
- **VI. ALGEBRAIC ITERATIVE METHODS** (103–113): A Asymptotic Source Extrapolation 103 ·
  B Chebyshev 104 · C Conjugate Gradient 107 · D Other Algebraic (Krylov: GMRES/BCGSTAB) 111
- **VII. ACCELERATION OF OTHER SCATTERING ITERATIONS** (113–126): A Upscattering 114 ·
  B Highly Forward-Peaked Scattering 118 · C Absorption-Emission in Radiative Transfer 121
- **VIII. ACCELERATION OF k-EIGENVALUE PROBLEMS** (126–137): A Power Iteration 127 ·
  B Shifted Power Iteration 129 · C Chebyshev for k 130 · D QD for k 131 · E Synthetic
  Acceleration for k 134 · F Summary 137
- **IX. DISCUSSION** (137–141) · Acknowledgements 141 · References 142–158 (≈ 302 entries)
  · Appendix A: List of Acronyms 158–159

**Open problems the Discussion enumerates** (pp. 139–140): (1) heterogeneity/aspect-ratio
degradation of consistent DSA, no known remedy; (2) parallel scalability of rapidly
convergent schemes (tied to scalable elliptic solves); (3) unstructured-grid acceleration;
(4) maturity of spatial multigrid and of QD-for-radiative-transfer; (5) the
CG/GMRES-outside-DSA trend (quoted in full at 2.7).

---

## Task 5 — DSA primaries reading map (for issue #2)

All four primaries are LOCAL in `scratch/literature/`. What the review says each
contributed:

- **Alcouffe (1977)** [63] — R.E. Alcouffe, "Diffusion Synthetic Acceleration Methods for
  the Diamond-Differenced Discrete-Ordinates Equations," *Nucl. Sci. Eng.* **64**, 344
  (preceded by the 1976 Trans. ANS note [57]). THE consistency breakthrough: "Alcouffe
  showed that if the discrete diffusion operator in the low-order problem is derived
  consistently with the high-order diamond-differenced S_N transport operator, then the
  resulting synthetic algorithm … is rapidly convergent for all spatial mesh thicknesses.
  This seminal result showed that it is possible to accelerate a discretized SI scheme
  using DSA in an efficient and robust manner" (pp. 16–17); "Prior to this work, it was
  not known whether an unconditionally fast acceleration method even existed for the
  discrete first-order transport equation" (p. 17). Also contains the nonlinear
  eigenvalue-ready variants (p. 58, p. 134) that became the DANTSYS production scheme,
  and the CMR/FMR/DSA head-to-head counts (CMR diverged/diverged/369; FMR 48/143/149;
  DSA 17/24/80 — p. 136). His source-correction form is algebraically equivalent to the
  four-step presentation the review teaches (p. 54).
- **Larsen (1982, Part I)** [83] — E.W. Larsen, "Unconditionally Stable Diffusion Synthetic
  Acceleration Methods for the Slab Geometry Discrete-Ordinates Equations. Part I:
  Theory," *Nucl. Sci. Eng.* **82**, 47. Generalizes Alcouffe's consistency principle
  into the **four-step procedure** applicable to any weighted-diamond/non-DD slab scheme
  (pp. 17, 58: "This same 'consistency' principle has been applied to other 1-D spatial
  differencing schemes by Larsen and McCoy [83, 85], whose method is sometimes called
  'Larsen's Four-Step' method"); the review's Eqs. (3.45)–(3.64) are its DD instance, and
  the GCB θ=3 application "has been reported [83]; … spectral radii are < c/3 for all
  values of c and cell thickness" (p. 73). The review also flags the deep detail that the
  four-step output depends on which algebraic form of the transport discretization you
  start from (which moments are taken) — p. 73 and the Figure 8 p₀-vs-p₁ moment
  sensitivity (pp. 76–77).
- **McCoy & Larsen (1982, Part II)** [85] — D.R. McCoy and E.W. Larsen, same title
  "…Part II: Numerical Results," *Nucl. Sci. Eng.* **82**, 64. The numerical-verification
  companion: demonstrates the Part-I theory across schemes; with [83] it anchors the
  claim that consistent-DSA "iterative convergence is extremely rapid: the iterative
  error reduces at least two orders of magnitude for every three DSA iterations" (p. 58).
- **Morel (1982)** [86] — J.E. Morel, "A Synthetic Acceleration Method for Discrete
  Ordinates Calculations with Highly Anisotropic Scattering," *Nucl. Sci. Eng.* **82**,
  34. The anisotropic-scattering extension: "Morel [86] showed that DSA applied to
  problems with highly anisotropic scattering will have a degraded performance, which can
  be partly remedied by accelerating both the flux and current using the solution of the
  low-order diffusion calculation, at least in 1-D planar geometry" (p. 17). Reads as the
  ancestor of the §VII.B story (per-moment eigenvalues ρ_n = Σ_sn/Σ_t; P₁-level
  acceleration of φ₀ AND φ₁), which the Manteuffel-Morel angular multigrid [164] later
  subsumed for the strongly forward-peaked regime.

Suggested reading order for #2 (slab DSA in ORPHEUS): review §§II.B + III.B (the tutorial
and the consistency mechanics) → Larsen Part I [83] (the general derivation to implement
against non-DD schemes, i.e. ORPHEUS's LD/DD family) → McCoy-Larsen Part II [85]
(verification targets: spectral-radius tables to reproduce as L1 evidence) → Alcouffe [63]
(DD specifics + eigenvalue variants) → Morel [86] (only when anisotropic acceleration
lands). For the GMRES-preconditioner path (#200), the operative results are Table 2
(p. 113) + the p. 140 design guidance; Ashby et al. [204] (*SIAM J. Numer. Anal.* **32**,
128 (1995)) is the citable journal version.

Supporting citations resolved from the bibliography (pp. 142–158): [10] Marchuk & Lebedev,
*Numerical Methods in the Theory of Neutron Transport*, 2nd rev. ed., Harwood (1986).
[25] Gol'din, *Zh. Vych. Mat. i Mat. Fiz.* **4**, 1078 (1964). [26] Kopp, NSE **17**, 65
(1963). [32] Gelbard & Hageman, NSE **37**, 288 (1969). [51] Reed, NSE **45**, 245 (1971).
[61] Lewis & Miller, Trans. ANS **23**, 202 (1976). [100] Derstine & Gelbard, Trans. ANS
**50**, 260 (1985). [110] Larsen & Miller, Trans. ANS **52**, 416 (1986). [139] Faber &
Manteuffel, in *Transport Theory, Invariant Imbedding, and Integral Equations*, Lect.
Notes Pure Appl. Math. 115, p. 37 (1989). [150] Cefus & Larsen, NSE **105**, 31 (1990).
[164] Morel & Manteuffel, NSE **107**, 330 (1991). [169] Adams & Martin, NSE **111**, 145
(1992). [180] Adams & Morel, NSE **115**, 253 (1993). [181] Adams & Wareing, Trans. ANS
**68**, 203 (1993). [190] Morel, Dendy & Wareing, NSE **115**, 304 (1993). [243] Ramoné,
Adams & Nowak, NSE **125**, 257 (1997). [257] Zika & Larsen, Trans. ANS **79**, 135
(1998). [259] Azmy, JCP **152**, 359 (1999). [260] Azmy, Wareing & Morel, Proc. M&C
Madrid, 1, 55 (1999). [279] Azmy, NSE **136**, 202 (2000). [298] Sanchez & Santandrea,
Proc. M&C Salt Lake City (2001). [300] Warsa, Wareing & Morel, Proc. M&C Salt Lake City
(2001). [287] Alcouffe, private communication (2001). Citation slip noted in the original:
§VIII.B cites SPI as "[131]" (p. 129, confirmed on the page image), but [131] is
Surzhikov's radiative-transfer preprint; the Wielandt reference in the bibliography is
[132] Sutton, NSE **98**, 169 (1988).

---

## Notation map and conflict flags (review → ORPHEUS)

- **c = Σ_s/Σ_t** — same as ORPHEUS scattering ratio. Heterogeneous predictor:
  min c ≤ ρ ≤ max c (3.31).
- **λ** — Fourier wave number **in mean-free-path units** (modes e^{iΣ_t λx}); their
  discrete **Λ = (2/(Σ_t h)) tan(Σ_t h λ/2)** (3.28). Flag: many later papers use
  unscaled λ (modes e^{iλx}); ORPHEUS docs should state the Σ_t-scaling explicitly when
  reproducing (2.50)/(3.65).
- **ω(λ)** = iteration eigenvalue (per-iteration mode attenuation); **σ = sup|ω|** =
  spectral radius = ORPHEUS ρ. β(λ) = SI eigenvalue when it appears inside DSA algebra
  (2.44); do not confuse with TSA's β (scattering-reduction parameter, Eq. (2.66)) or the
  CG β step-scalar (6.30) — three unrelated β's in one paper.
- **Iterate superscripts**: ℓ+1/2 = post-sweep, ℓ+1 = post-acceleration (2.22)–(2.34);
  k-eigenvalue chapter uses n for outers, s for spectrum/group sub-iterations (§VIII.D).
- **Operators — the big conflict flag**: the review's **A ≡ I − L⁻¹S** (Eq. (1.27),
  angular-flux level) and **T ≡ I − X_SI** (Eq. (6.4), scalar-flux level, X_SI =
  K₀L⁻¹(Σ_s/2)) are both *sweep-preconditioned* fixed-point operators. **Neither is the
  ORPHEUS honest matvec A = L + C − S − B** (un-preconditioned, typed-operator algebra).
  Mapping: ORPHEUS's GMRES-on-honest-matvec-with-sweep-preconditioner solves exactly the
  review's Tφ = η̂ system (Eqs. (6.5), (6.21)–(6.26)); the review's L = "leakage plus
  collision" (1.4) ≈ ORPHEUS L + C (their L absorbs the collision term; ORPHEUS splits
  streaming L from collision C and carries the boundary operator B explicitly).
- **Their P = I + MS** (1.32) is the preconditioner acting on the Richardson residual;
  DSA's M = ED⁻¹K (1.40) with K = P₁ projection, D = low-order operator, E = prolongation
  — a clean match to ORPHEUS frame/projection vocabulary (K, E are the analysis/
  reconstruction faces of an angular frame; the four-step D is "apply K to the discrete
  L − S", which is the whole consistency doctrine in one sentence).
- **Slab μ** = ORPHEUS `mu_x`. **Cylinder (p. 80)**: their **ξ** = axial-direction cosine
  = ORPHEUS `mu_z`; their **μ** = "the ANGLE between the outward radial direction and the
  radial-plane projection of the direction of particle flight" — an *angle, not a
  cosine* (flag: most cylinder literature, incl. ORPHEUS, uses the cosine pair; ORPHEUS
  `mu_x` = √(1−ξ²)·cos(their μ)). **Sphere (p. 78)**: their μ = cos θ against the outward
  radial direction = ORPHEUS `mu_x` convention for 1-D sphere.
- **Their "unstable/divergent"** is used interchangeably for σ > 1 (p. 54: "divergent or
  unstable (σ_DSA > 1)"); the docs should keep ORPHEUS's sharper usage but note the
  review's equivalence.
- **k-eigenvalue**: A = L⁻¹P (8.5) — here L is the full loss operator (streaming +
  collision − scattering), i.e. ORPHEUS A = L+C−S with B folded into BCs; their P =
  production = ORPHEUS F. Their dominance ratio r = ORPHEUS's DR; SPI shift κ obeys
  (8.25): 0 < κ−k₁ ≪ |κ−k₂|, with the warning that shifted inners are "nearly-critical
  [… very highly scattering] problems" so "the fixed-source acceleration techniques …
  become even more important" (p. 129), and in multigroup the shifted fission source
  "introduces an effective upscattering" (p. 130).

---

*Survey complete. Deliverable for corpus Phase I; do not merge into
`documentation_corpus_architecture.md` — the main agent owns that edit. Extraction
working file: `scratch/literature/adams_larsen_2002_paged.txt` (page-marked pdftotext
dump; journal p. N marker `=== p.N ===`).*
