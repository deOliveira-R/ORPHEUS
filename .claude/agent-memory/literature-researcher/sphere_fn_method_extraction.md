---
name: Sphere F_N method extraction (Siewert-Thomas 1986)
description: Sphere F_N derivation for Branch 1 — verdict: Siewert-Thomas 1986 covers slabs+spheres in one unified formalism with anti-symmetric BC; 1G is clean special case (drop matrix structure to scalars). No additional paper required.
type: project
---

# Verdict

**Siewert-Thomas 1986 (NSE 94, 264-270, DOI 10.13182/NSE86-A17269) is sufficient.**
The paper does slabs AND spheres simultaneously in a unified F_N formalism.
The 1G special case is a clean reduction (collapse all 2x2 matrix structure
to scalars). **No additional paper needs to be acquired.** All four supporting
references (Siewert-Benoist Part I, Grandjean-Siewert Part II, KLL 1974,
Mitsis 1963) are locally available in `scratch/literature/` for cross-checks.

# Sphere = slab with one sign flip — the load-bearing fact

The sphere problem on r∈[0,R] is solved by the substitution Ψ(r,μ)·r → ψ(x,μ)
on x∈[-R,R] with **anti-symmetric** parity ψ(x,μ) = -ψ(-x,-μ), versus the slab's
**symmetric** parity ψ(x,μ) = ψ(-x,-μ). In Siewert-Thomas this is one BC change:

- **Slab BC** (Eq. 4): Ψ(-a,μ) = Ψ(a,-μ), μ∈[0,1]
- **Sphere BC** (Eq. 46, p. 268): Ψ(-a,μ) = **-**Ψ(a,μ), μ∈[0,1]

Quoted from p. 268: *"the critical sphere problem requires only that Eq. (4) be
changed to read [Eq. 46], and that we interpret a as the critical radius, we
have incorporated the relevant minus sign in our developed equations, and thus
we list in Table IV our results for critical radii for spheres."*

This means **every F_N matrix entry in the slab derivation reuses verbatim
for sphere except the relative sign of the two exp(-2a/ν) terms**. In Eq.
(43)/(46) the pattern `[X_α + Y_α exp(-2a/ν_β)] a_α = 0` becomes
`[X_α - Y_α exp(-2a/ν_β)] a_α = 0` for sphere (or equivalently, the
right-hand normalization vector flips sign). Mechanically the simplest
encoding is: pass a `geometry_sign ∈ {+1, -1}` (slab=+1, sphere=-1) into
the F_N matrix assembler.

# Sphere F_N matrix entries — the equations

Verbatim from Siewert-Thomas, applied with geometry_sign s ∈ {+1, -1}:

**Critical condition** (Eq. 39 generalized): for ν_β a discrete or collocation
point in (0, 1/σ) ∪ (0, 1):

  Σ_{α=0}^N [X_α(ν_β) + s · Y_α(ν_β) exp(-2a/ν_β)] a_α = 0    (39a)

**X_α and Y_α scalar entries** (1G reduction, σ=1, C=c, c_11=c, others=0;
det CL terms collapse; σ_2_f drops):

  Y_α(ν) = c · F_α(ν)                                       (32a, 1G)
  X_α(ν) = -c · F_α(-ν)                                     (36a, 1G)

where (Eqs. A.2-A.4):

  F_α(ξ) = ∫_0^1 μ · P_α(2μ-1) · dμ / (μ + ξ)               (A.2)
  G_α(ξ) = P · ∫_0^1 μ · P_α(2μ-1) · dμ / (ξ - μ)           (A.3)

with G_α(ξ) = -F_α(-ξ) for ξ∈(0,1) (symmetry relation, p. 269).

**Forward recursion for F_α (Eqs. A.5/A.6)**, real ξ ∉ [-1.001, -1]:

  F_0(ξ) = 1 - ξ ln(1 + 1/ξ)                                (A.5)
  F_{α+1}(ξ) = -[(2α+1)(2ξ+1) F_α(ξ) − α F_{α-1}(ξ) + δ_{α,0} + δ_{α,1}] / (α+1)   (A.6)

For real ξ ∉ [-1.001, 0.001]: backward recursion (Eqs. A.7-A.10) — keeps
F_α stable for ξ in the F_N collocation grid which avoids both [-1,0] and
neighborhoods of zero.

**Collocation points** (Eq. 38a, 1G case where σ=1, ℵ=1):

  ξ_β = (1/2)[1 + cos(βπ/(N+1))],   β = 1, 2, ..., N        (38a)

i.e. shifted Chebyshev-of-first-kind nodes on (0,1), plus β=0 contributes
the dispersion-root equation (39a evaluated at ν=ν_0).

**Critical condition (1G sphere)**: det M(a) = 0, where M is the (N+1)×(N+1)
matrix whose β-th row (β=0,...,N) has entries

  M_{β,α} = X_α(ξ_β) + s · Y_α(ξ_β) · exp(-2a/ξ_β)

with s = -1 for sphere, ξ_0 ≡ ν_0 the discrete eigenvalue (root of Λ(ν)=0,
Eq. 18 of Westfall-Metcalf or Eq. 7 of Grandjean-Siewert: Λ(z) = 1 + (cz/2)
∫_{-1}^{1} dμ/(μ-z)). Newton's method on det M(a) = 0 gives the critical
radius a = R_c.

# Critical condition: clean form

For 1G isotropic sphere with c > 1:

  **det[X(ξ_β) - Y(ξ_β) exp(-2R/ξ_β)] = 0**   (sphere)
  **det[X(ξ_β) + Y(ξ_β) exp(-2R/ξ_β)] = 0**   (slab)

The slab and sphere F_N differ by exactly **one global sign on the
exp-coupled block**. This is the architecture point: a single Branch 1
SymPy module can serve both geometries via one parameter.

# X-function reusability

**Yes — slab and sphere share the same Wiener-Hopf X-function.** The X(z)
defined in Siewert-Benoist Part I Eq. (18),

  X(z) = (1/(1-z)) · exp[(1/π) ∫_0^1 arg Λ⁺(τ) · dτ/(τ-z)]

depends only on the dispersion function Λ(z) = 1 + (cz/2) ∫_{-1}^{1} dμ/(μ-z),
which is a property of the medium (parameter c), not the geometry. Both slab
and sphere F_N use the **same** X(z) for the residual / B-coefficient definitions
(though for the F_N approximation itself, X(z) doesn't appear explicitly — it's
embedded in the "exact" reference values used to validate F_N convergence).
Reusable: any X-function code already in `core/` for slab F_N transfers verbatim.

# Flux reconstruction

For the sphere with anti-symmetric ψ on x∈[-a,a], the boundary expansion
(Eq. 25, modified) becomes:

  ψ(a, μ) = Σ_{α=0}^N a_α P_α(2μ-1),  μ ∈ [0,1]              (25')

The interior reconstruction follows from the integral form of the
Boltzmann equation (Case-Zweifel §3.5 / Mitsis 1963):

  ρ(r) = (1/r) · ψ(x=r, μ) integrated against full-range eigenfunctions

Concrete formula: from Siewert-Thomas's slab Eq. (40)/(41) with
A_{α,2} → -a_α (sign flip for sphere) and the standard a_0 = -1
normalization (or any chosen power level), the sphere flux is
reconstructed by inverse-transforming the half-range expansion.
**Defer flux reconstruction in implementation** — k_eff is the
primary deliverable and tests against Sood `Ua-1-0-SP` reference
need only the critical radius, not flux shape.

# 1G special case — clean or hard?

**Clean.** The 2G→1G reduction is mechanical:

| 2G object              | 1G reduction                                |
| ---------------------- | ------------------------------------------- |
| σ = diag{σ, 1}, σ>1    | σ = 1 (scalar)                              |
| C = 2×2 transfer matrix| C = c (scalar)                              |
| Θ(μ) = diag{θ(μ),1}    | Θ(μ) = 1                                    |
| Λ(z) = I + zC ∫dμ Θ/(μ-z) | Λ(z) = 1 - cz·atanh(1/z)                |
| Ω(z) = CΛ(z)C⁻¹        | Ω(z) = Λ(z)                                 |
| ℵ = 1 or 2 (eigenval count) | ℵ = 1 always (single ν_0)              |
| det CL(ν_β)            | det collapses; "det" → identity             |
| det ω(ν) = ω(ν)        | scalar                                      |

The X_α / Y_α formulas in Eqs. (32)-(36) all collapse cleanly: matrix
products become scalar products, the 2-component vectors a_α become
scalar a_α, and the system size (Eqs. 42a-c with three scalar equations)
collapses to **one** scalar equation per collocation point (just Eq. 42a).
Final 1G F_N matrix is (N+1)×(N+1).

# Implementation feasibility per stage

| Stage                              | Feasibility | Notes                                               |
| ---------------------------------- | ----------- | --------------------------------------------------- |
| F_α(ξ), G_α(ξ) closed form         | SymPy clean | `F_0 = 1 - ξ·log(1 + 1/ξ)`, recursion symbolic     |
| F_α recursion (Eqs. A.5/A.6)       | SymPy clean | Forward stable for ξ > 0 (Bessel-style stability)  |
| Dispersion function Λ(z)           | SymPy clean | `1 - c·z·atanh(1/z)`                                |
| Discrete eigenvalue ν_0            | SymPy/mpmath| Newton root-find: `Λ(ν_0) = 0`, ν_0 > 1 for c<1, imaginary for c>1 — use mpmath for c>1 |
| Collocation grid ξ_β (Chebyshev)   | SymPy clean | `cos(βπ/(N+1))`                                     |
| F_N matrix assembly                | SymPy clean | Pure algebra in F_α(ξ_β)                            |
| Critical condition det M(a) = 0    | mpmath      | 1D root-find on R; matrix entries contain `exp(-2a/ξ)` — closed form in SymPy, numeric Newton in mpmath |
| Flux reconstruction (deferred)     | -           | Defer; not needed for k_eff cross-check             |

**SymPy walls: none.** All special functions are elementary
(`log`, `atanh`, `cos`, `exp`). The F_α recursion is polynomial-rational
in ξ and stable.

**mpmath fallback needed only for**:
1. Imaginary discrete eigenvalue ν_0 when c > 1 (SymPy's `nsolve` works
   but mpmath is faster for high precision).
2. The 1D root-find for R_c at user-requested precision (e.g. 10 sig-figs
   to match Sood's published 2.4248249802).
3. Matrix-determinant evaluation at large N where symbolic det is slow —
   mpmath's `mpf` matrix det (LU-based) is robust.

**Special-function obstacles: none.** Unlike cylinder F_N (Bickley-Naylor
Ki_n series, modified Bessel K_0/K_1, Westfall-Metcalf 1972 machinery),
sphere F_N stays in the slab function basis (shifted Legendre on [0,1]
with logarithmic kernels). This is the **easiest** geometry to extend.

# Notation map (Siewert-Thomas → ORPHEUS Branch 1)

| S-T 1986          | ORPHEUS Branch 1 1G slab/sphere              |
| ----------------- | -------------------------------------------- |
| a                 | half-thickness (slab) / radius R (sphere)    |
| μ                 | direction cosine, ∈[-1,1]                    |
| Ψ(x,μ)            | angular flux; for sphere: ψ(r,μ) = Ψ(x,μ)/r  |
| σ                 | =1 in 1G                                     |
| C, c_ij           | = c (scalar mean secondaries)                |
| ν_0               | discrete eigenvalue (Case spectrum)          |
| Λ(ν), Ω(z)        | dispersion function (scalar in 1G)           |
| F_α(ξ)            | reused verbatim — Eq. A.2                    |
| ξ_β (collocation) | Chebyshev nodes Eq. 38a                      |
| s ∈ {+1, -1}      | geometry_sign: slab=+1, sphere=-1            |

# Action items

**No additional paper needs to be acquired.** All references are locally
available:

- Primary: `scratch/literature/Siewert-Thomas(1986)On Two-Group Critical Problems in Neutron-Transport Theory.pdf`
- Slab F_N theory baseline: `Siewert-Benoist(1979)The FN Method ... Part I_ Theory and Applications.pdf`
- Slab F_N numerics: `Grandjean-Siewert(1979)... Part II_ Applications and Numerical Results.pdf`
- Reference 1G sphere benchmark: `Sood Foster Parsons (1999)Analytical Benchmark Test Set for Criticality Code Verification.pdf` → `Ua-1-0-SP`, R_c = 2.4248249802 mfp at c=1.30
- 1G sphere derivation (alternative — Mitsis-style for cross-check): `Mitsis(1963) Transport Solutions to the Monoenergetic Critical Problems.pdf`
- 1G sphere benchmark machinery (alternative — Wiener-Hopf Fredholm, KLL): `Kaper-Lindeman-Leaf(1974)Benchmark Values for the Slab and Sphere Criticality Problem in One-Group Neutron Transport Th.pdf`

**Recommendation for method-implementer**: build the Branch 1 SymPy
sphere F_N as a parameterization of the slab F_N module via
`geometry_sign ∈ {+1, -1}`. The only sphere-specific code is one
sign in the matrix assembly. F_α recursion, X-function machinery,
collocation grid, dispersion solver — all reuse from slab F_N
verbatim. Validate against Sood `Ua-1-0-SP` (2.4248249802 mfp,
c=1.30) at N=5 (expect ≤4-sig-fig agreement per Grandjean-Siewert
Table V conclusion) and N=10 or N=30 for full 6-sig-fig match
against Siewert-Thomas Table IV reference column.
