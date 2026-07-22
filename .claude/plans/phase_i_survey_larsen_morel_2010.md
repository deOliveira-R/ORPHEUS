# Phase-I survey — Larsen & Morel (2010), "Advances in Discrete-Ordinates Methodology"

**Full citation (CrossRef-verified 2026-07-21):** Larsen, E.W., and Morel, J.E.,
"Advances in Discrete-Ordinates Methodology," Chapter 1 in Y. Azmy and
E. Sartori (eds.), *Nuclear Computational Science: A Century in Review*,
Springer, Dordrecht (2010), **pp. 1–84**. DOI `10.1007/978-90-481-3411-3_1`
(CrossRef: published online 2009-12-24; container "Nuclear Computational
Science"; pages "1-84").

**Source examined:** the LOCAL file
`/Users/rodrigo/git/nuclear/ORPHEUS/scratch/literature/Nuclear Computational Science - A Century in Review.pdf`
(whole book, 476 PDF pages). **Page map: book page N = PDF page N+12**
(chapter 1 = PDF 13–96). Equation-bearing pages were read visually
(pdftotext mangles the Springer math encoding); prose/references were
cross-checked against the text extraction. All page numbers below are
**book pages**. Zotero was not queried this session (tools not exposed;
the brief named the local file directly — no annotations checked).

---

## Task 1 — §3.2 table row

| Column | Verdict | Evidence |
|---|---|---|
| **S_N page count** | **84 pp. total** (confirmed): main text pp. 1–75, Authors' Note + 126 refs pp. 76–82, author bios pp. 83–84. CrossRef records "1-84". | Chapter footer p. 1; CrossRef. |
| **Cylinder?** | **Essentially NO.** No cylindrical S_N equations are ever written. Curvilinear treatment = **1-D spherical geometry only**: §1.2 pp. 7–9 (Eqs. 1.20–1.26, the original discretization) + §1.5.1 pp. 36–37 (the flux-dip story, Eqs. 1.92–1.93) ≈ **5 pages of spherical, 0 pages of cylindrical**. 2-D (r,z) cylindrical geometry gets ONE sentence: "This scheme has been generalized to 2-D cylindrical geometry [38]" (p. 37, [38] = Morel & Montry 1984). The only other cylinder mention is an application aside — "pencil-beam sources in 2-D cylindrical and 3-D geometries" for charged particles (p. 20). Grep: `cylind` hits on exactly 2 pages (20, 37). | pp. 7–9, 20, 36–37. |
| **Adjoint of S_N?** | **NO — never constructed.** The adjoint of the discretization (à la Bell & Glasstone §6.2b) does not appear. "Adjoint" occurs only as: (a) adjoint *solutions* used as a Monte Carlo variance-reduction input in the oil-well-logging application (pp. 21–23: "the second is to use deterministic methods to provide adjoint calculations for automatic variance reduction in Monte Carlo calculations [118]"; ATTILA adjoint → MCNP weight windows, "a factor of 1,000 speedup," p. 22); (b) "self-adjoint" as an operator-symmetry PROPERTY — self-adjoint angular flux codes for charged particles (p. 49, [112]) and the "self-adjoint and positive-definite (SAPD)" property of the integral-moments operator for CG (p. 71, [107]). No adjoint equations, no duality identity, no discrete-adjoint commutation discussion anywhere. | pp. 21–23, 49, 71. |
| **Eigenvalue posed independently?** | **NO — the k-eigenvalue problem is never posed.** The LBE (Eq. 1.1, p. 2) carries the fission production term inline (χ/4π ∫∫νΣ_f ψ) as part of a **fixed-source, time-dependent** posing; no k division, no criticality eigenproblem, no power iteration. The word "criticality" never appears. The nearest approaches: §1.8.3 accelerates "the fission source in steady-state subcritical or time-dependent super-critical neutronics calculations [73]" (p. 60) — still fixed-source with fission; and the k-eigenvalue exists ONLY as bibliography entry [126] (Warsa, Wareing, Morel, McGhee, Lehoucq 2004, "Krylov subspace iterations for deterministic k-eigenvalue calculations," NSE 147:26), cited once inside the Krylov-adoption citation list on p. 63. | pp. 2, 60, 63, 82. |
| **Verification?** | **ZERO.** Grep across the whole chapter: "verification"/"verify"/"verified" = **0 hits**; "manufactured" = **0**; "validat*" = **0**. "Benchmark" = **1 hit** (p. 22) and it is an application anecdote, not code verification: "Monte Carlo calculations for both tools were performed with the MCNP code [72] to provide benchmark results" (ATTILA vs MCNP for logging-tool detector responses). "Test problems" = 1 hit (p. 37, the Walters–Morel [60] angular-scheme accuracy comparison). No analytical benchmarks, no MMS, no convergence-order measurement methodology, no discussion of how any of the review's accuracy claims would be *evidenced* by a reader. | Whole-chapter greps; pp. 22, 37. |

**Headline for the corpus thesis:** even the *modern* review monograph —
explicitly written to cover what the textbook canon missed — contains **no
verification content at all**, no adjoint-of-the-discretization, no
independent eigenvalue posing, and no cylindrical S_N formulation. The
mini-book niche survives contact with this reference intact, and gains: the
review is instead the **best single source for M4b convergence theory**
(Task 2), which the three textbooks also lack.

---

## Task 2 — M4b harvest: convergence theory

This is where the review earns its place. §1.2 (pp. 11–13), §1.4.4
(pp. 30–36), §1.8 (pp. 54–63), §1.9 (pp. 63–74) together constitute the
only textbook-form treatment of S_N iterative + asymptotic convergence
theory in the surveyed canon.

### 2.1 Spectral radius of Source Iteration

Definition (p. 12), two forms — note (1.33b) is the *computable* estimator:

$$\rho = \lim_{\ell\to\infty} \frac{\lVert \psi - \psi^{(\ell)}\rVert}{\lVert \psi - \psi^{(\ell-1)}\rVert} \tag{1.33a}$$

$$\rho = \lim_{\ell\to\infty} \frac{\lVert \psi^{(\ell)} - \psi^{(\ell-1)}\rVert}{\lVert \psi^{(\ell-1)} - \psi^{(\ell-2)}\rVert} \tag{1.33b}$$

The headline result (p. 12):

$$\rho = \frac{\Sigma_s}{\Sigma_t} = c \equiv \text{scattering ratio} \tag{1.34}$$

> "Therefore, the SI scheme applied to one-group problems always converges;
> but if c is arbitrarily close to unity (the probability of absorption is
> arbitrarily small) and the system becomes increasingly thick, SI will
> converge arbitrarily slowly. (As c → 1, increasingly many particles will
> undergo increasingly many scattering events during their histories.)" (p. 12)

Physical interpretation of the iterate (p. 11): with zero initial guess,
"the SI angular flux estimate after ℓ sweeps is, physically, the angular
flux due to particles that have experienced at most ℓ−1 scattering events."

Multigroup: SI always converges, but high upscatter probability → ρ near
unity (p. 12). Thermal-radiation outer iteration (§1.3.1, p. 17): inner
ρ = c = Σ*_s/Σ*_t; outer

$$\rho_0 = \upsilon \int_0^\infty \chi\!\left(\frac{\Sigma_a^*}{\Sigma_a^* + \tau}\right)\mathrm{d}E, \qquad \tau = \frac{1}{c\,\Delta t^k}\ \ (1.44\text{b}) \tag{1.46}$$

with the warning "ρ₀ ≈ 0.9999" is possible (p. 17). (Here the linearized
radiative-transfer equation (1.43) "has the form of the steady-state
neutron transport equation, with Σ*_a acting as the fission cross section,
υ acting as the number of neutrons per fission, χ acting as the fission
spectrum" — p. 16.)

### 2.2 Fourier analysis (§1.8.1, pp. 55–57) — the full machinery

Model: infinite homogeneous medium, no discretization, Eqs. (1.165)
(transport) and (1.166) (SI). Exact error equations (1.167)–(1.168) for
δφ^(ℓ), δψ^(ℓ−1/2). Ansatz (p. 55):

$$\delta\phi^{(\ell-1)}(x) = \int_{-\infty}^{\infty} a^{(\ell-1)}(\lambda)\, e^{i\Sigma_t\lambda x}\,\mathrm{d}\lambda \tag{1.169a}$$

$$\delta\psi^{(\ell-1/2)}(x,\mu) = \int_{-\infty}^{\infty} b^{(\ell-1/2)}(\lambda,\mu)\, e^{i\Sigma_t\lambda x}\,\mathrm{d}\lambda \tag{1.169b}$$

which yields per-mode algebra

$$(i\lambda\mu + 1)\, b^{(\ell-1/2)}(\lambda,\mu) = \frac{c}{2}\, a^{(\ell-1)}(\lambda) \tag{1.170a}$$

$$b^{(\ell-1/2)}(\lambda,\mu) = \frac{c}{2}\frac{1}{1+i\lambda\mu}\, a^{(\ell-1)}(\lambda) \tag{1.171}$$

and the **iteration eigenvalue** (p. 56):

$$\omega(\lambda) = \frac{c}{2}\int_{-1}^{1} \frac{\mathrm{d}\mu'}{1+i\lambda\mu'} = \frac{c}{\lambda}\tan^{-1}(\lambda) \tag{1.172b}$$

$$\lVert\delta\phi^{(\ell)}\rVert \approx \rho^\ell A \ \ (1.174), \qquad
\rho = \sup_{-\infty<\lambda<\infty} |\omega(\lambda)| \tag{1.175}$$

$$\rho = \sup_{-\infty<\lambda<\infty}\left(\frac{c}{\lambda}\tan^{-1}\lambda\right) = c \tag{1.176}$$

"which is attained for λ ≈ 0" (p. 57). The **shape of the slow modes** is
the DSA motivation: for λ ≈ 0,
b^(ℓ+1/2)(λ,µ) ≈ (c/2)(1 − iλµ)ω^ℓ(λ)a^(0)(λ), so

> "the most slowly converging modes depend nearly linearly on µ." (p. 57)

Scope and reliability of the tool (p. 57) — quote worth carrying whole:

> "The Fourier analysis can be extended to multidimensional problems, to
> fully discrete (in space, angle, and energy) S_N problems for infinite
> homogeneous (or spatially periodic) media on Cartesian grids, and to
> iteration strategies beyond Source Iteration [110]. This analysis makes
> it possible to predict – with relative ease and great accuracy – the rate
> of convergence of transport iteration schemes, before implementing them
> in test codes. It also provides a theoretical foundation for iteration
> schemes, making it possible to understand how a scheme works (if it
> works), or why it fails (if it fails)."

And: "For the SI scheme, the Fourier analysis predicts that the spectral
radius ρ = c, even for fully discrete S_N codes. The accuracy of this
prediction has been observed in many calculations [110]." (p. 57)
([110] = Adams & Larsen 2002, PNE 40:3 — the designated "recent
comprehensive review" of iterative methods, p. 54.)

### 2.3 DSA and its consistency theory (§1.8.2, pp. 57–59)

Algorithm: sweep (1.177) → exact error transport equation (1.178) → replace
by its diffusion approximation

$$-\frac{\mathrm{d}}{\mathrm{d}x}\frac{1}{3\Sigma_t}\frac{\mathrm{d}}{\mathrm{d}x}\,\delta\phi^{(\ell-1/2)}(x) + \Sigma_a\,\delta\phi^{(\ell-1/2)}(x) = \Sigma_s\left[\phi^{(\ell-1/2)}(x) - \phi^{(\ell-1)}(x)\right] \tag{1.179}$$

$$\phi^{(\ell)}(x) = \phi^{(\ell-1/2)}(x) + \delta\phi^{(\ell-1/2)}(x) \tag{1.180}$$

Spectral radius (p. 58): "The DSA spectral radius for this case is
**approximately 0.23c**, which is bounded less than unity for all
0 ≤ c ≤ 1." (The classic 0.2247c; the review rounds.) Mode-complementarity
explanation (p. 58):

> "The transport sweep strongly attenuates scalar flux errors that vary
> rapidly in space, i.e., high-frequency errors. However, this sweep
> attenuates low-frequency errors that vary slowly in space by only a
> factor of c. … In contrast, the diffusion step almost completely
> attenuates low-frequency errors (because such errors have an angular
> shape that is primarily linear in µ) but does basically nothing to
> high-frequency errors. … The minimum level of attenuation occurs for
> intermediate-frequency errors."

**Consistency theorem history** (p. 59) — the load-bearing paragraph:

> "The concept of DSA existed long before 1968 [6], but the synthetic
> method for discrete problems was originally seen to be unstable for
> problems in which the cell thickness exceeded roughly a mean free path
> [11]. Alcouffe [23] made the DSA method practical for the
> diamond-differenced S_N equations by showing that if the spatial
> discretization of the diffusion equation was chosen to be *consistent*
> with the spatial discretization of the S_N equations, the instability was
> eliminated, and one obtained an algorithm that was unconditionally
> effective. Soon afterward, Alcouffe's ideas were extended and generalized
> to nondiamond schemes [34, 35, 99, 108, 110]."

([6] Kopp 1963; [11] Gelbard & Hageman 1969; [23] Alcouffe 1977 NSE 64:344;
[34]/[35] Larsen / McCoy & Larsen 1982 — the Fourier-analyzed
unconditionally-stable slab DSA pair.)

Formal definitions (p. 59): "an acceleration scheme is said to be
'**unconditionally effective**' when the accelerated spectral radius is
bounded less than unity. Also, an acceleration scheme is said to be
'**unconditionally efficient**' when the computational execution time
associated with the accelerated scheme is always much smaller than that of
the unaccelerated scheme." Effective ≠ efficient: fully consistent
diffusion discretizations of advanced (DFEM) S_N schemes are nonstandard
and expensive → "an effective-but-inefficient DSA algorithm [116]"; modern
codes use "partially-consistent" discretizations [64] or solve the fully
consistent system approximately [61, 67], accepting "a degradation of the
accelerated spectral radius" for efficiency (p. 59).

Quasidiffusion caveat (p. 59): Gol'din's method [8] is "rapidly convergent
but is not strictly an acceleration method, because it produces converged
solutions that differ from the discrete S_N solution due to an extra
truncation error."

Anisotropic scattering (p. 59): standard DSA leaves higher Legendre moments
unaccelerated → ineffective for forward-peaked problems; accelerating the
current (n = 1) helps in 1-D [36] but "In multidimensional calculations,
acceleration of the currents becomes unstable as the anisotropy of the
scattering increases [68]. For 1-D calculations, an angular multigrid
method has been developed, which is unconditionally effective and efficient
with anisotropic scattering [62]. However, this approach has had only
limited success in multidimensional calculations [93]."

**The multidimensional deficiency** (§1.8.4, pp. 62–63):

> "For many years, Alcouffe's consistent DSA approach appeared to yield an
> unconditionally effective DSA algorithm. However, Azmy [87] has recently
> shown that the DSA method becomes increasingly ineffective in
> heterogeneous multidimensional calculations as the jumps in cross
> sections across material discontinuities increase in magnitude. It is
> likely that all DSA-like methods will exhibit this same deficiency. …
> Fortunately, this deficiency has been overcome by recasting DSA as a
> preconditioner and using it in conjunction with preconditioned Krylov
> methods."

Krylov-side statement of the same fact (p. 71): "Only one eigenvalue of the
iteration matrix very near unity can ruin the performance of an
acceleration scheme, but a single anomalously large eigenvalue associated
with the system being solved will generally have little effect on a Krylov
method." Also p. 63: "an inconsistently discretized DSA method, when recast
as a preconditioner for a Krylov method, can produce an unconditionally
effective and efficient S_N solution technique [63]."

### 2.4 Asymptotic thick-diffusion-limit theory (§1.4.4, pp. 30–36)

The Larsen–Morel–Miller lineage, presented as the derivation itself
([49] = LMM 1987 JCP 69:283; [52] = Larsen & Morel 1989 JCP 83:212
+ corrigendum JCP 91:246 (1990); [66] = Larsen 1992 NSE 112:336 review;
[86, 88, 105] extend it).

Motivation — why truncation-order analysis is the wrong instrument (p. 30):
"an nth-order scheme would satisfy ‖ψ_exact − ψ_Δx‖ = O(τⁿ), τ ≪ 1 …
In slab geometry, the DD and SC schemes are second-order, while the LD
method is third-order." But: "a truncation error analysis does not provide
useful information for these types of problems. … To determine if accurate
results can be obtained with a mesh that is optically thick but resolves
the spatial scale length of the solution, it is necessary to perform a
discrete asymptotic analysis [49, 113]." Also the multi-D order pathology
(p. 30): "in multidimensional geometries the S_N solutions have singular
characteristics, across which the solution is not smooth [20] … computer
experiments have shown [33] that because of the singular characteristics,
the order of convergence of the DD scheme in x,y-geometry depends on the
definition of the error norm."

The admissibility statement (p. 31):

> "A spatial discretization of the S_N equations is of practical use for
> diffusive problems if it possesses the *optically thick diffusion limit*
> [49, 52, 66, 86, 88, 105]. Such a discretization scheme will yield
> accurate results for diffusive problems if the spatial mesh cells are
> thin with respect to a diffusion length … even if these cells are thick
> with respect to a mean free path."

The scaling (Eqs. 1.77–1.78): Σ_t → Σ_t/ε, Σ_a → εΣ_a, Q → εQ (diffusion
equation invariant; S_N equations not). Scaled S_N (1.79). Result:

$$\psi_n(x) = \frac{\phi(x)}{2} + O(\varepsilon) \tag{1.80}$$

via the expansion ψ_n = Σ_i ε^i ψ_n^(i) (1.81) and the hierarchy (1.82):

- i = 0 (1.83): leading order isotropic, ψ_n^(0) = φ^(0)/2 (1.85) —
  **requires** Σ w_n = 2 (1.84);
- i = 1 (1.86): ψ_n^(1) = φ^(1)/2 − (µ_n/2Σ_t) dφ^(0)/dx (1.88) —
  **requires** Σ µ_n w_n = 0 (1.87);
- i = 2 (1.89): has a solution only under the **solvability condition**
  (1.91), which IS the diffusion equation
  0 = d/dx (1/3Σ_t) dφ^(0)/dx − Σ_a φ^(0) + Q — **requires**
  Σ µ_n² w_n = 2/3 (1.90).

The quadrature moment conditions (1.84)/(1.87)/(1.90) are the *angular
admissibility* hypotheses of the spatial limit theorem — worth surfacing
as such in the mini-book.

Thick diffusion limit defined (p. 34): Σ_t → O(1/ε) with ψ_n → φ/2; "the
total cross section becomes unbounded, yet the S_N solution limits to an
O(1) diffusion solution." The discrete question: "*What happens if this
same asymptotic analysis is applied to the spatially discretized S_N
equations?*" — with the subtlety that Σ_tΔx → ∞ as ε → 0, so naive
accuracy reasoning gives the wrong answer either way (p. 34).

**Scheme verdicts** (p. 34) — the table no textbook has:

> "Some schemes are accurate in the thick diffusion limit; others are not.
> For example, the Step-Characteristic (SC) scheme fails as ε → 0 (the SC
> solution → 0). The Diamond-Difference (DD) scheme fails unless all the
> diffusive regions of the problem have isotropic incident boundary fluxes
> (in the presence of nonisotropic boundary fluxes, DD solutions become
> corrupted by unphysical spatial oscillations). LD-like schemes perform
> successfully in the thick diffusion limit in 1-D geometries. LD methods
> also perform well in multi-D geometries with triangular (2-D) or
> tetrahedral (3-D) spatial grids, but they fail in quadrilateral (2-D) or
> hexahedral (3-D) grids. (However, bilinear-discontinuous methods work
> well for quadrilateral grids and trilinear-discontinuous methods work
> well for hexahedral grids.)"

(This matches the ORPHEUS memory record `multi_d_ld_closure.md`: simplex-P1
fails on quads/hexes per Adams 2001 [105]; the 2^d-corner bilinear/
trilinear family is the fix.)

**Unresolved boundary layers** (p. 35) — an open-problem statement:

> "the conclusions so far are that no known differencing scheme is
> completely adequate to model unresolved boundary layers accurately. For
> example, LD methods are generally inaccurate in the first cell
> (containing the boundary layer) within the thick diffusive region, and
> they incorrectly predict that the flux exiting the diffusive region is
> isotropic. Generally, to be certain that a discrete solution is accurate,
> all spatial boundary layers must be adequately resolved by the spatial
> grid."

**Fokker–Planck limit** (p. 35): charged particles have a different
asymptotic limit — Σ_t = O(ε⁻¹) with mean scattering cosine
µ̄₀ = 1 − O(ε); "Space-angle discretization schemes have also been
successfully analyzed in this asymptotic limit [113]" ([113] = Pautz &
Adams 2002). "Ensuring that the discretized S_N equations limit to a valid
discretization of the Fokker–Planck equation is primarily related to the
treatment of anisotropic scattering rather than the spatial differencing
scheme." (Contrast: the diffusion limit is spatial-differencing-dominated.)

### 2.5 Krylov convergence theory (§1.9, pp. 63–74)

- Minimum-polynomial argument that x = A⁻¹b ∈ K_d(A,b) (Eqs. 1.197–1.200,
  p. 64); Krylov-vector degeneracy: "these vectors become less linearly
  independent with increasing m … A^j b approaches the fundamental
  eigenvector" (p. 65) → orthogonalization; weighting-space taxonomy:
  Ritz–Galerkin W_m = K_m (CG, minimizes A-norm of error, Eqs. 1.206–1.207)
  vs least-squares W_m = AK_m (GMRES/MINRES, minimizes L₂ residual,
  Eqs. 1.208–1.210); restart caveat: "convergence of the restarted process
  is not guaranteed for general matrices. However, the restarted GMRES
  method *is* guaranteed to converge if A is positive-definite, i.e., if
  x^T A x > 0 for all x ≠ 0" (pp. 66–67).
- Condition number κ = σ_max/σ_min (1.211, p. 67) governs CG; "The
  condition number is not particularly relevant to the convergence of the
  generalized minimum-residual (GMRES) method. More relevant factors are
  the distribution of the eigenvalues of A in the complex plane, and …
  the condition number of the matrix of eigenvectors" (p. 67). Clustering:
  "If any single property of the coefficient matrix is desirable with a
  Krylov method, it is that its eigenvalues should be clustered away from
  zero." (p. 67)
- Preconditioner theory for transport (p. 68): "the best preconditioners
  move the smallest eigenvalues significantly away from zero, while
  leaving the largest eigenvalues relatively unaffected [124]." The
  anti-pattern (p. 69): preconditioners that also spread the largest
  eigenvalues [121, 122] are "analogous to unstable acceleration schemes
  that strongly attenuate the low-frequency Fourier error modes but weakly
  amplify the high-frequency Fourier error modes. We emphasize that not
  every acceleration scheme that attenuates low-frequency … can be
  expected to be effective when recast as a preconditioned Krylov method."
- **Why sweep-preconditioning works** (p. 70, on Eq. 1.218
  (I − ½L⁻¹Σ_sP)ψ = L⁻¹Q with L = µ∂_x + Σ_t (1.216), P = ∫dµ (1.217)):
  "the analytic operator on the left side of Eq. (1.218) represents the
  integral transport operator, which is bounded, whereas the differential
  transport operator is unbounded. … Furthermore, the integral transport
  operator is a compact perturbation of the identity operator, hence many
  (but not all) of the eigenvalues of the discretized S_N version of the
  integral operator will be clustered about unity." Positive-definiteness
  of the discretized (1.218) "has been observed (but not rigorously proven
  [124])" → restarted-GMRES guarantee (p. 70).
- Peierls form (1.219) (I − ½PL⁻¹Σ_s)φ = PL⁻¹Q: smaller unknown vector;
  SAPD-able [107]; "The SAPD property has been observed to be preserved
  for Eq. (1.219) on orthogonal meshes with a wide variety of S_N spatial
  discretization schemes [82], but it is not preserved with the
  linear-discontinuous spatial discretization scheme on unstructured
  tetrahedral meshes [124]." (p. 71)
- DSA-preconditioned moments system (1.220)
  (I + D⁻¹Σ_s)(I − ½PL⁻¹Σ_s)φ = (I + D⁻¹Σ_s)PL⁻¹Q, D = −∂_x(1/3Σ_t)∂_x + Σ_a
  (1.221): "the analog of diffusion-synthetic acceleration, and it
  (presumably) results in rapid convergence under almost all practical
  conditions [124]. The deficiency of DSA in multidimensions with large
  discontinuities in the cross sections is strongly diminished when it is
  recast as a preconditioned Krylov method [124]." (p. 71). pp. 72–74 give
  the operational sweep-level recipe (source-vector build + operator
  action) proving "only sweeps, diffusion solves, and vector additions and
  subtractions are required" (p. 74).

### 2.6 Admissibility criteria inventory (feeds the first-class "admissibility" section)

Collected criteria the review states as requirements on a scheme, with locus:

1. **Positivity** (pp. 4–7, 10): implicit time differencing "is guaranteed
   to yield a positive solution. (This is not the case for other time
   discretizations.)" (p. 5); multigroup and Cartesian S_N are "inherently
   positive"; "The positivity of the S_N approximation for 1-D planar and
   other Cartesian geometries does not automatically apply in curvilinear
   geometries" (p. 7); DD positive only for sufficiently small cells in
   1-D, "never guaranteed" in multi-D (p. 10). "A desirable feature of a
   discretization for the LBE is that the resulting approximate solution
   should be positive – or nearly so." (p. 4)
2. **WD spatial weights** (p. 10): "the constants α_{n,j} are constrained
   by the goals of accuracy, stability, and symmetry about µ = 0."
3. **Angular-coefficient constant-solution condition** (p. 8): β_{n+1/2}
   "is defined so that Eq. (1.23a) admits a constant solution for the
   infinite medium configuration in the presence of a constant source."
4. **Quadrature moment conditions** for the diffusion limit: Σw = 2 (1.84),
   Σµw = 0 (1.87), Σµ²w = 2/3 (1.90).
5. **Thick-diffusion-limit possession** (p. 31, verdicts p. 34) — the
   scheme-level admission test for diffusive problems.
6. **Thermal-radiation scheme requirements** (p. 17): "spatial differencing
   schemes must be highly damped in strongly absorbing media, accurate in
   optically thin regions, have the thick diffusion limit …, and behave
   well with unresolved spatial boundary layers."
7. **DSA consistency** (p. 59): low-order operator discretization consistent
   with the S_N discretization ⇒ unconditional effectiveness (Alcouffe).
8. **Rebalance grid admissibility** (p. 13): fine-mesh unstable; coarse
   enough stable; too coarse ineffective [56].
9. **Galerkin quadrature invertibility** (p. 41): "D and M should be
   inverses of one another" (Gauss: M = D⁻¹ automatic; otherwise generate
   moment-dependent weights (1.100) or invert M directly).
10. **Companion/Galerkin weight positivity** (p. 75): "These quadrature
    weights should be positive because (among other things) positivity
    indicates a stable spherical-harmonic interpolation. … a
    multidimensional Galerkin quadrature set beyond S₂ with positive
    weights for all moments has never been constructed."
11. **SAPD preservation** (p. 71): discretization should preserve the
    analytic SAPD property of the integral-moments operator to unlock CG.
12. **Angular-differencing linear-preservation** (p. 37): the diamond-in-
    angle scheme "does not preserve solutions that are linear in µ" —
    the M-M weighted diamond (1.93) restores consistency with the cosine
    location (see Task 5a).

---

## Task 3 — M2 harvest: symptom → cause → fix triples

| # | Symptom | Cause | Fix | Locus |
|---|---|---|---|---|
| 1 | Negative fluxes + "unphysical oscillatory behavior" in DD solutions | 1-D: cell widths too large; multi-D: "positive diamond-difference solutions are never guaranteed; oscillatory behavior is always possible, for any spatial grid" | Weighted-diamond α_{n,j} ≠ 0 (worked "in more strongly absorbing reactor shields"); nonlinear negative-flux fixup (worked in "diffusive and weakly absorbing reactor cores"); or LD/DFEM (p. 27: "much more accurate and robust … than finite-difference methods") | p. 10; p. 27 |
| 2 | Fixup-induced solver trouble: "the resulting equations can be difficult to solve even when nonlinear solution techniques are used" | Set-negative-outflows-to-zero fixup [71] is "highly nonlinear and nondifferentiable"; bites when "a problem contains both highly absorptive and highly scattering regions" | Emerging "repair" paradigm — do not fix up during iteration; repair after convergence [123] (Shashkov & Wendroff 2004) | p. 74 |
| 3 | **Discrete-ordinates flux dip** — "an erroneous suppression in the flux at the center of the sphere" (recognized early 1960s, not eliminated until early 1980s) | THREE compounding causes (pp. 36–37): (i) starting-direction equation discretized in the curvilinear-like form (1.92) — "thought to make the discretization … consistent with that of the other directions. However, it actually contributed to truncation errors that enhanced the flux dip" [22]; (ii) specular-reflection center BC — "even though the angular flux at the center of a sphere is rigorously isotropic … This incorrectly allowed the angular flux at r = 0 to be anisotropic"; (iii) "The diamond-in-angle equation is inconsistent with the location of the quadrature cosines within each angular bin. As a result, the diamond-in-angle scheme does not preserve solutions that are linear in µ" | (i) discretize the slab-geometry form of the starting-direction equation (1.25); (ii) set ALL center angular fluxes equal to the starting-direction value [26]; (iii) Morel–Montry angular weighted diamond (1.93) "consistent with the location of the cosine in each angular bin" [38]. "When all three of these measures were combined, the resulting angular discretization scheme eliminated the flux dip [38]." | pp. 36–37 |
| 4 | SI converges arbitrarily slowly | c → 1 with optically thick system (diffusive); or high upscatter probability (multigroup) | DSA (§1.8.2); upscatter two-grid acceleration [69]; fission-source acceleration [73]; LMFG for radiative transfer [42, 50] | pp. 12, 57–62 |
| 5 | Synthetic (diffusion) acceleration **diverges** | Cell thickness > ~1 mfp with an inconsistent low-order discretization [11] | Alcouffe consistency [23]; or recast as Krylov preconditioner (immune: "a single anomalously large eigenvalue … will generally have little effect on a Krylov method," p. 71) | pp. 59, 63, 71 |
| 6 | Fine-mesh rebalance **diverges**; coarse-mesh rebalance ineffective | "For most problems, fine-mesh rebalance is unstable (divergence occurs), while coarse-mesh rebalance with a sufficiently coarse grid is stable and provides acceleration. However, if the coarse grid is too coarse, then coarse-mesh rebalance becomes ineffective [56]" | Choose the coarse grid carefully (costly setup); historically superseded by DSA | p. 13 |
| 7 | DD solutions "corrupted by unphysical spatial oscillations" in thick diffusive regions | Nonisotropic incident boundary fluxes on a diffusive region (DD lacks the thick diffusion limit in that case) | LD-like schemes (1-D, simplex multi-D); bilinear/trilinear discontinuous on quads/hexes | p. 34 |
| 8 | SC solution → 0 (vanishes) in diffusive regions | SC "fails as ε → 0" — no thick diffusion limit | Use LD-family / DFEM schemes there | pp. 25, 34 |
| 9 | Wrong first-cell answer + spuriously isotropic exiting flux at a diffusive interface | Unresolved boundary layer (few-mfp transition volume not meshed) | No known scheme fix — "all spatial boundary layers must be adequately resolved by the spatial grid" | p. 35 |
| 10 | DSA ineffective for forward-peaked (charged-particle) problems | Higher-order Legendre flux moments unaccelerated and slow to converge | Accelerate current too (1-D) [36]; multi-D current acceleration goes UNSTABLE with increasing anisotropy [68]; 1-D angular multigrid unconditionally effective+efficient [62], "only limited success" in multi-D [93] | pp. 20, 59 |
| 11 | Multi-D DSA effectiveness degrades in heterogeneous problems | Cross-section jumps across material discontinuities (Azmy [87]); "likely that all DSA-like methods will exhibit this same deficiency" | Recast DSA as a Krylov preconditioner [124] | pp. 62–63, 71 |
| 12 | Krylov iterates stagnate / degrade numerically | Krylov vectors lose linear independence (A^j b → fundamental eigenvector); restarted GMRES not guaranteed for general matrices | Orthogonalize (Arnoldi); ensure positive-definiteness (then restarted GMRES is guaranteed); CG three-term recursion when SPD | pp. 65–67 |
| 13 | Negative discrete scattering sources despite positive discrete fluxes | Truncated Legendre/polynomial interpolation of the angular flux "can be negative at some points, even though the discrete values themselves are positive" | With Gauss quadrature the negativities are provably SMALL ("polynomial interpolation at the Gauss points is known to be stable"); for non-Gauss sets, Galerkin/moment-dependent weights (1.100) — but "no guarantee that the weights generated in this way will be positive" | pp. 39–40 |
| 14 | **Ray effects** (the review's "most significant of all S_N deficiencies") | Discrete-direction angular collocation; "the exact solutions to such problems usually have an extremely anisotropic angular dependence" | "Perhaps surprisingly, very little has been accomplished during the past 40 years" (p. 36) [21, 24]; P_N no cure ("they rarely provide the correct solutions to problems for which S_N solutions exhibit ray effects"); angular FEM coupling "does not appear promising" [57, 119]; "the accurate elimination of ray effects will probably require angular adaptation" [85, 92, 111, 117] | pp. 36, 75 |
| 15 | Accuracy LOSS when replacing the starting-direction flux by a discontinuous angular treatment | "the starting-direction flux is computed (by Eq. (1.25)) with greater accuracy than the other directions; hence, significant accuracy is actually lost if the starting-direction flux plays no role in the angular derivative treatment" | Keep the starting-direction sweep; if using discontinuous-in-angle, use quadratic-continuous in the FIRST angular cell + LD in the rest [60] | p. 37 |
| 16 | Numerical diffusion in energy (charged-particle spectra smeared) | Multigroup-like (diamond) treatment of continuous slowing down | LD/DFEM in energy: "much less numerically diffusive than the step method and much less oscillatory than the diamond method" [47] | pp. 20, 51 |

---

## Task 4 — Section structure + the review's own scope statement

### Scope/gap claim (p. 1, §1.1) — verbatim

> "In 1968, Bengt Carlson and Kaye Lathrop published a comprehensive review
> on the state of the art in discrete-ordinates (S_N) calculations [10]. …
> Since 1968, several books and reviews on general numerical methods for
> S_N simulations have been published [32, 46, 71], but **none of these
> covers the advanced work done during the past 20 years**."

([10] Carlson & Lathrop 1968; [32] Sanchez & McCormick 1982 NSE 80:481;
[46] Marchuk & Lebedev 1986; [71] Lewis & Miller 1984.) And the design
principle (p. 1): "The specific purpose of this chapter is to describe how
the field of S_N calculations has matured **through the lens of three
important physical problems** that can be simulated today but could not be
realistically simulated in 1968." — i.e., the review is *application-pulled*
(radiative transfer, charged particles, well logging), which explains the
reactor-side gaps (no k-eigenvalue, no cylinder): its driver problems are
not reactor problems.

### Full outline with page ranges and equation spans

| Section | Pages | Equations | Content |
|---|---|---|---|
| 1.1 Introduction | 1–2 | — | 1968 baseline; gap claim; roadmap |
| 1.2 Basic Concepts | 2–14 | 1.1–1.36 | LBE (1.1–1.4); planar (1.5–1.8); implicit time (1.9–1.12); multigroup (1.13–1.17); S_N slab (1.18–1.19); **1-D spherical** (1.20–1.26); spatial FD + DD/WD (1.27–1.30); SI (1.31–1.32); spectral radius (1.33–1.34); rebalance (1.35–1.36) |
| 1.3 Three Challenging Physical Problems | 14–23 | 1.37–1.50 | — |
| 1.3.1 Thermal Radiation Transport in the Stellar Regime | 14–18 | 1.37–1.49 | Nonlinear T-coupling; Newton linearization → effective-fission structure (1.43–1.44); inner/outer ρ (1.46); equilibrium-diffusion limit (1.47–1.49) |
| 1.3.2 Charged-Particle Transport | 18–21 | 1.50 | Boltzmann–Fokker–Planck splitting; forward-peaked issues |
| 1.3.3 Oil-Well Logging Tool Design | 21–23 | — | ATTILA vs MCNP anecdote; adjoint-for-VR |
| 1.4 Advances in Spatial Discretizations | 23–36 | 1.51–1.91 | — |
| 1.4.1 Characteristic Methods | 23–25 | 1.52–1.55 | SC derivation; 1-D SC ≡ a WD scheme |
| 1.4.2 Linear Discontinuous Method | 25–27 | 1.56–1.61 | LD balance + slope equations; lumped LD; corner balance |
| 1.4.3 Nodal Methods | 27–29 | 1.62–1.73 | Transverse integration; CCN/LLN |
| 1.4.4 Solution Accuracy in the Thick Diffusion Limit | 30–36 | 1.74–1.91 | The asymptotic derivation + scheme verdicts + boundary layers + FP limit |
| 1.5 Advances in Angular Discretizations | 36–45 | 1.92–1.124 | — |
| 1.5.1 Angular Derivatives | 36–37 | 1.92–1.93 | Flux dip: 3 causes, 3 fixes; M-M weighted diamond |
| 1.5.2 Anisotropic Scattering | 38–45 | 1.94–1.124 | Legendre closure result; Galerkin quadrature; M/D matrices; delta/straight-ahead invariance; extended transport correction |
| 1.6 Advances in Fokker–Planck Discretizations | 45–52 | 1.125–1.157 | — |
| 1.6.1 The Continuous-Scattering Operator | 45–49 | 1.125–1.147 | Angular FP operator; moment-preserving discretizations; effective moments Σ_em = Σ_r,tr[(N−1)N − m(m+1)]/2 (1.147) |
| 1.6.2 The Continuous-Slowing-Down Operator | 49–52 | 1.148–1.157 | CSD sweeps in energy; LD-in-energy vs diamond |
| 1.7 Advances in Time Discretizations | 52–54 | 1.158–1.164 | DFEM-in-time; nodal-in-time |
| 1.8 Advances in Iteration Acceleration | 54–63 | 1.165–1.193 | — |
| 1.8.1 Fourier Analysis | 55–57 | 1.165–1.176 | ρ_SI = c; slow modes linear in µ |
| 1.8.2 Diffusion-Synthetic Acceleration | 57–59 | 1.177–1.180 | ρ_DSA ≈ 0.23c; consistency; effectiveness vs efficiency |
| 1.8.3 DSA-Like Methods for Outer Iteration Acceleration | 60–62 | 1.181–1.193 | LMFG; fission-source & upscatter acceleration; rank-one vs full-rank operators |
| 1.8.4 A Deficiency in Multidimensional DSA and DSA-Like Methods | 62–63 | — | Azmy heterogeneity result; Krylov rescue |
| 1.9 Krylov Methods | 63–74 | 1.194–1.222 | — |
| 1.9.1 The Central Theme of Krylov Methods | 63–67 | 1.194–1.210 | Krylov spaces; CG/GMRES taxonomy |
| 1.9.2 Convergence and Preconditioning of Krylov Methods | 67–69 | 1.211–1.214 | κ; clustering; preconditioner theory |
| 1.9.3 Applying Krylov Methods to the S_N Equations | 69–74 | 1.215–1.222 | First-order vs integral vs Peierls vs DSA-preconditioned forms; sweep-level recipe |
| 1.10 Future Challenges | 74–75 | — | DFEM on nonorthogonal grids across 3 limits; positivity/repair; unstructured parallel sweeps; curved faces; optimal preconditioned Krylov; pencil beams; positive Galerkin weights; ray effects/angular adaptation |
| Authors' Note + References | 76–82 | — | 126 refs, chronological bands (1900–1964, 1965–1969, …, 2000–2005) |
| Author bios | 83–84 | — | Larsen; Morel |

**Structural-precedent observations for the mini-book:** (a) the
organization is *by discretized variable* (space → angle → energy[FP] →
time) then *by solver layer* (stationary acceleration → Krylov), with a
motivating-applications chapter up front and an open-problems chapter at
the end; (b) "Basic Concepts" replays the whole 1968-era pipeline
(time → energy → angle → space → iteration) in 12 pages so every later
section can diff against it; (c) each advance is presented WITH its failure
history (what was wrong before, who fixed it, what remains broken) — the
closest existing analogue of the corpus's symptom→cause ambition, but
scattered as narrative rather than tabulated; (d) references are grouped
chronologically, not topically — useful as a dated map of the field,
useless for lookup-by-topic.

---

## Task 5 — Bonus sweeps

### 5a. Curvilinear angular-coefficient conventions (feeds the M5 crosswalk)

**L&M-2010 names the angular redistribution coefficients β, not α.**
The conservative-form spherical S_N balance (p. 8), after multiplying by r²:

$$\mu_n \frac{\partial}{\partial r} r^2 \psi_n(r) + \frac{r}{w_n}\left[\beta_{n+1/2}\psi_{n+1/2}(r) - \beta_{n-1/2}\psi_{n-1/2}(r)\right] + \Sigma_t(r) r^2 \psi_n(r) = \sum_{n'=1}^{N} \Sigma_s(r,\mu_n,\mu_{n'})\, r^2 \psi_{n'}(r) w_{n'} + \frac{r^2 Q(r)}{2} \tag{1.23a}$$

$$\beta_{n+1/2} = -2\sum_{n'=1}^{n} \mu_{n'} w_{n'} \tag{1.23b}$$

with diamond-in-angle ψ_n = ½(ψ_{n−1/2} + ψ_{n+1/2}) (1.24) and the
starting-direction equation (µ = −1, ∂ψ/∂µ bounded):

$$-\frac{\partial \psi_{1/2}(r)}{\partial r} + \Sigma_t(r)\,\psi_{1/2}(r) = \sum_{n'=1}^{N} \Sigma_s(r,-1,\mu_{n'})\,\psi_{n'}(r) w_{n'} + \frac{Q(r)}{2} \tag{1.25}$$

Convention analysis (extends the normalization table in agent memory
`curvilinear_sweep_directness_ruling.md`):

- (1.23b) ⇔ recursion **β_{n+1/2} = β_{n−1/2} − 2µ_n w_n, β_{1/2} = 0**,
  and β_{N+1/2} = 0 automatically by quadrature symmetry. This is
  **IDENTICAL to the Bailey–Morel–Chang sphere α-recursion (BMC Eq. 11)**
  — same ×2, same Σw = 2 normalization, same divisor structure
  (r/w_n) — the symbol is simply renamed β.
- The **minus sign is absorbed into the definition** (partial sums of µw
  over the inward hemisphere are negative, so β ≥ 0) and the
  redistribution term enters with a **+** sign in (1.23a). Flag per the
  standing notation rule: "the sign before the redistribution depends on
  whether α absorbs the minus sign."
- ⚠ **Symbol-collision hazard**: L&M-2010 uses **α_{n,j} for the SPATIAL
  weighted-diamond parameter** (Eq. 1.30, p. 10):
  ψ_{n,j} = ((1+α_{n,j})/2)ψ_{n,j+1/2} + ((1−α_{n,j})/2)ψ_{n,j−1/2},
  α = 0 ⇒ diamond. In ORPHEUS/Bailey usage α is the ANGULAR coefficient.
  Crosswalk: L&M α_{n,j} ↔ ORPHEUS spatial WD weight τ via τ = (1+α)/2
  (for µ_n > 0); L&M β_{n+1/2} ↔ ORPHEUS/BMC α_{m+1/2}.
- Design rationale stated for β (p. 8): "defined so that Eq. (1.23a)
  admits a constant solution for the infinite medium configuration in the
  presence of a constant source" — the constant-solution-preservation
  admissibility criterion.
- The **Morel–Montry angular weighted diamond** (p. 37) is transcribed:

$$\psi_n(r) = \left(\frac{\mu_n - \mu_{n-1/2}}{w_n}\right)\psi_{n+1/2}(r) + \left(\frac{\mu_{n+1/2} - \mu_n}{w_n}\right)\psi_{n-1/2}(r) \tag{1.93}$$

  ≡ BMC Eq. 15 with the M-M weight τ_m = (µ_m − µ_{m−1/2})/w_m (BMC
  Eq. 42) ≡ Bailey Eq. 74 ≡ the ORPHEUS implementation. L&M-2010 provides
  the *narrative* justification ("consistent with the location of the
  cosine in each angular bin", restores linear-in-µ preservation) but NOT
  the asymptotic derivation.
- µ mapping: L&M-2010 spherical µ = Ω·r/r (Eq. 1.20) = ORPHEUS `mu_x`
  (streaming/radial cosine, spherical convention per L-001).

### 5b. Angular-differencing asymptotics (Bailey–Morel–Chang lineage)

**The review contains NO asymptotic diffusion-limit analysis of ANGULAR
differencing.** §1.4.4's asymptotic theory is applied to SPATIAL schemes
only; §1.5.1's angular-derivative treatment is the flux-dip narrative
(causes + M-M fix + Lathrop-2000 comparison [100] + Walters-Morel [60]),
with no ε-scaling of the angular closure. Bailey, Morel & Chang,
"Asymptotic Diffusion-Limit Accuracy of S_N Angular Differencing Schemes"
(NSE 165:149, 2010) is **absent from the 126-entry bibliography** (chapter
finalized ~2009; BMC published the same year). Consequence for the corpus:
the review documents the spatial half of the diffusion-limit consistency
story (LMM 1987 [49] / L&M 1989 [52] / Larsen 1992 [66] / Adams 2001
[105]); the angular half exists ONLY in the BMC journal paper (local:
`scratch/literature/Bailey-Morel-Chang(2010)…pdf`; extraction:
agent memory `space_angle_discretization_separability.md`). A mini-book
section that states BOTH per-axis conditions (spatial ⊗ angular,
factorized but jointly required) would exceed every existing review
including this one. The only space-angle joint asymptotics in the review
is the Fokker–Planck limit sentence (p. 35, [113] Pautz & Adams 2002).

Adjacent angular-axis content worth reusing: the Galerkin-quadrature
operator algebra (pp. 40–45) — S⃗ = **MΣD**ψ⃗ (1.102) with
D_{m,n} = P_m(µ_n)w_n (1.104), M_{n,m} = ((2m+1)/2)P_m(µ_n) (1.106),
Gauss ⇒ **M = D⁻¹** (1.107) ⇒ S⃗ = D⁻¹ΣDψ⃗ (1.108): "the S_N scattering
matrix represents a similarity transformation of the Legendre scattering
matrix" (p. 41) — i.e., the S_N/P_{N−1} slab equivalence in operator form.
ORPHEUS mapping: L&M **D** ≡ HarmonicFrame *analysis* face (ψ → moments),
L&M **M** ≡ *reconstruction* face (moments → ψ); the M = D⁻¹ requirement is
the frame-invertibility (perfect-reconstruction) condition. Also the
straight-ahead invariance (1.115)–(1.124): Σ*₀ − MΣ*D = Σ₀ − MΣD (1.124),
so the extended transport correction (α = Σ_N, traditional near-optimal
choice, p. 45) "leaves the S_N solution invariant" while cutting c — the
convergence/solution-invariance separation.

### 5c. Bibliography entries most relevant to DSA (feeds issue #2)

Chronological, as numbered in the chapter (pp. 76–82). Core lineage:

- [6] Kopp HJ (1963) "Synthetic method solution of the transport equation." NSE 17:65. — pre-history.
- [8] Gol'din VYa (1964) quasi-diffusion. USSR Comp Math Math Phys 4:136. — nonlinear cousin (converged answer differs).
- [11] Gelbard EM, Hageman LA (1969) "The synthetic method as applied to the S_N equations." NSE 37:288. — the instability discovery.
- [16] Reed WH (1971) "The effectiveness of acceleration techniques for iterative methods in transport theory." NSE 45:245. — thin-grid-only performance.
- [23] **Alcouffe RE (1977)** "Diffusion synthetic acceleration methods for the diamond-differenced discrete-ordinates equations." NSE 64:344. — THE consistency paper. (Local PDF: `scratch/literature/Alcouffe(1977)….pdf`.)
- [34] **Larsen EW (1982)** "Unconditionally stable DSA methods for the slab geometry discrete-ordinates equations. Part I: theory." NSE 82:47; [35] McCoy DR, Larsen EW (1982) Part II: numerical results. NSE 82:64. (Local PDFs: `Larsen(1982)… Part II.pdf` + `Larsen(1982)….pdf`.)
- [36] Morel JE (1982) "A synthetic acceleration method for discrete ordinates calculations with highly anisotropic scattering." NSE 82:34. (Local: `Morel(1982)….pdf`.)
- [42] Morel JE, Larsen EW, Matzen MK (1985) JQSRT 34:243 — radiative-diffusion synthetic acceleration (LMFG root); [50] Larsen EW (1988) "A grey transport acceleration method…" JCP 78:459.
- [44] Khalil H (1985) "A nodal diffusion technique for synthetic acceleration of nodal S_N calculations." NSE 90:263.
- [56] Cefus GR, Larsen EW (1990) "Stability analysis of coarse mesh rebalance." NSE 105:31.
- [61] Wareing T, Larsen EW, Adams ML (1991) DFEM DSA slab/x-y (M&C Pittsburgh); [67] Wareing T (1992) LA-12425-T.
- [62] Morel JE, Manteuffel TA (1991) "An angular multigrid acceleration technique for the S_N equations with highly forward-peaked scattering." NSE 107:330; [93] Pautz, Morel, Adams (1999) multi-D angular multigrid (M&C Madrid); [68] Adams ML, Wareing TA (1993) "DSA given anisotropic scattering, general quadratures, and multi-dimensions." TANS 68:203.
- [63] Ashby SF, Brown PN, Dorr MR, Hindmarsh AC (1991) "Preconditioned iterative methods for discretized transport equations" (M&C Pittsburgh). — inconsistent-DSA-as-preconditioner works.
- [64] **Adams ML, Martin WR (1992)** "Diffusion synthetic acceleration of discontinuous finite element transport iterations." NSE 111:145. — partially-consistent DSA.
- [69] Adams BT, Morel JE (1993) two-grid upscatter acceleration. NSE 115:253; [84] Adams BT, Morel JE (1997) fission+upscatter (M&C Saratoga Springs).
- [70] Morel JE, Dendy JE Jr, Wareing TA (1993) "Diffusion accelerated solution of the 2-D S_N equations with bilinear-discontinuous differencing." NSE 115:304.
- [73] Morel JE, McGhee JM (1994) "A fission-source acceleration technique for time-dependent even-parity S_N calculations." NSE 116:73.
- [74] Brown PN (1995) "A linear algebraic development of DSA for three-dimensional transport equations." SIAM J Numer Anal 32:179.
- [82] Ramone GL, Adams ML, Nowak PF (1997) "A transport synthetic acceleration method…" NSE 125:257; [103] Zika MR, Adams ML (2000) TSA with opposing reflecting BCs. NSE 134:159; [121]/[122] Hanshaw et al. (2003) stretched/filtered TSA.
- [87] **Azmy YY (2002)** "Unconditionally stable and robust adjacent-cell diffusive preconditioning of weighted-difference particle transport methods is impossible." JCP 182:213. ⚠ printed under the "1995–1999" band — misfiled. Companion positive results: [94] Azmy (1999) JCP 152:359; [99] Azmy (2000) NSE 136:202 (adjacent-cell preconditioners).
- [98] Wareing, Morel, McGhee (1999) DSA for DFEM space+time (M&C Madrid).
- [108] Sanchez R, Chetaine A (2001) "A synthetic acceleration for a two-dimensional characteristic method in unstructured meshes." NSE 136:122.
- [110] **Adams ML, Larsen EW (2002)** "Fast iterative methods for discrete-ordinates particle transport calculations." PNE 40:3. — the designated comprehensive iterative-methods review (local: `scratch/literature/Adams-Larsen(2002)….pdf`).
- [116] **Warsa JS, Wareing TA, Morel JE (2002)** "Fully consistent DSA of linear discontinuous S_N transport discretizations on unstructured tetrahedral meshes." NSE 141:236. — effective-but-inefficient.
- [124] **Warsa JS, Wareing TA, Morel JE (2004)** "Krylov iterative methods and the degraded effectiveness of DSA for multidimensional S_N calculations in problems with material discontinuities." NSE 147:218. — the modern synthesis (deficiency + Krylov rescue).
- [126] Warsa JS, Wareing TA, Morel JE, McGhee JM, Lehoucq RB (2004) "Krylov subspace iterations for deterministic k-eigenvalue calculations." NSE 147:26.

Reference-list hygiene notes: [87] misfiled by band (above); [88] Adams,
Wareing, Walters printed "(1988)" but NSE 130:18 is **1998**; [51] Baker &
Koch printed "(1988)" but NSE 128:312 is **1998** — three dating slips a
citing document should silently correct.

### 5d. Incidental finds worth keeping

- **SC-as-WD equivalence** (p. 25): "For 1-D problems, the SC method can be
  formulated as a weighted-diamond scheme, i.e., of the form described by
  Eqs. (1.28) and (1.30)" — one WD parameterization spans DD/SC/step in
  1-D (supports a single ORPHEUS `CellUpdate` family).
- **Truncation orders on record** (p. 30): slab DD and SC are 2nd-order,
  LD 3rd-order; multi-D order is norm-dependent because of singular
  characteristics ([20] Kellogg 1974, [33] Larsen 1982 NSE 80:710) — a
  ready-made caveat for any MMS order-verification chapter.
- **Legendre-truncation invariance** (p. 38): a degree-L cross-section
  expansion with a degree-L flux gives the SAME scattering source as the
  exact kernel — "In this case, the convergence of the cross-section
  expansion is irrelevant. This powerful result is not widely appreciated."
- **Future-challenges list** (pp. 74–75) doubles as the review's own gap
  table: DFEM accuracy+robustness across the thin/thick-absorptive/
  thick-diffusive trio on nonorthogonal meshes; positive+monotone
  second-order schemes without nondifferentiable fixup (repair paradigm
  [123]); parallel unstructured sweeps [104, 114]; curved cell faces [109];
  optimal preconditioned Krylov; pencil beams [125]; positive
  multidimensional Galerkin weights; ray-effect elimination via angular
  adaptation.
