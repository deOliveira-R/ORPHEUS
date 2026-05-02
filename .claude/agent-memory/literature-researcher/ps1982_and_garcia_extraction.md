---
name: ps1982_and_garcia_extraction
description: Pomraning-Siewert 1982 and Garcia 2017/2019/2021 PDF extractions for Plan 2 Variant α verification (vacuum sphere precursor + multi-region stable-P_N). Records exactly what each paper provides as L1 reference and — critically — what each paper does NOT provide for the k-eigenvalue Plan 2 use case.
type: project
---

# Pomraning-Siewert 1982 + Garcia 2017/2019/2021 — extraction memo

PDFs read in full from `/workspaces/ORPHEUS/scratch/literature/`. For
Plan 2 Variant α verification (continuous-µ multi-bounce sphere kernel
with vacuum or specular BC, eventually multi-region with reflection).

---

## TL;DR — by paper

**PS-1982 (4-pp Note, JQSRT 28).** Derives the integral form of the
**fixed-source** transfer equation for a homogeneous sphere with
isotropic scattering, allowing simultaneous specular reflection
coefficient α, diffuse reflection coefficient β·χ(µ), and external
illumination K(µ). Eq. (21) is a self-contained Fredholm-2 integral
equation for `r·I(r)` whose kernel is `[E_1(|r−x|) − E_1(r+x)]` plus
two BC-correction kernels `αF_1(r,x) + βF_2(r,x)` that vanish in the
vacuum limit α=β=K=0. **No k-eigenvalue results, no critical-radius
calculations, no numerical solutions** — this is a derivation Note
only. The mathematical path is **not** the cosh-even-extension path
of Sanchez 1986; it integrates Eq. (8a)+(8b) over µ and adds (the
"two halves of the sphere" trick Mitsis used). **This is structurally
independent of Sanchez 1986** at the level of derivation but the
final kernel is identical to Sanchez 1986 ω₁=0 limit.

**Garcia 2021 (J. Comp. Phys. 433, 19 pp).** Builds an iterative
region-by-region stable-P_N solver for the multi-region sphere with
internal sources and externally incident angular flux. **Subcritical
only — `c_k < 1`, criticality stated as future work.** Solves three
test cases on a 3-region sphere (R=3/5/7 cm) at orders N=1,3,9,19,59,199.
Reports converged scalar flux φ_k(r) and current J_k(r) tables to
5 significant figures. **No k-eigenvalue values reported in this paper.**
Garcia's machinery is structurally independent of Variant α — spherical
harmonics expansion + Marshak projection vs integral operator along
characteristics — but only the **fixed-source benchmark** is cross-checkable;
the k-eigenvalue cross-check Plan B3 originally targeted is **not
available from this paper**.

**Garcia 2017 (JQSRT 196).** A 4-pp Note giving the *particular* P_N
solution for a position-and-angle-dependent source — the building
block for handling a non-trivial `S(r,µ)` in Eq. (1). Cited as
Ref [21] in Garcia 2021. Pure formula derivation, no numerical
results.

**Garcia 2019 (JCP 405).** Numerically stable P_N solution for a
spherical **shell** (single annular region). Reformulates the angular
redistribution term as a heterogeneous source so the streaming-only
homogeneous P_N is plane-geometry-like. Particular solution requires
solving a system of N+1 Volterra integral equations of the second kind.
Cited as Ref [15] in Garcia 2021 — the "shell solution" reused for
intermediate and outer regions of the multi-region method.

**Garcia 2019b (JCTT 47).** Stable P_N for the **exterior** of a sphere
(transmission-type problem). Establishes that ill-conditioning in shell
problems comes from modified spherical Bessel functions of the **third**
kind (k_n), whereas a pure sphere uses only first-kind (i_n). Ref [38]
in Garcia 2021.

---

## Part 1 — PS-1982 (priority for Plan A2)

### Citation
Pomraning, G.C. and Siewert, C.E., "On the integral form of the
equation of transfer for a homogeneous sphere," J. Quant. Spectrosc.
Radiat. Transfer **28**(6), 503–506 (1982).
DOI `10.1016/0022-4073(82)90016-4`. Received 25 January 1982.

### Problem statement (Eq. 1, Eq. 2)

Equation of transfer in optical units, r ∈ (0, R], µ ∈ [−1, 1]:

$$
\left\{\mu \frac{\partial}{\partial r} + \frac{1}{r}(1-\mu^2)\frac{\partial}{\partial\mu} + 1\right\} I(r,\mu)
   = \frac{\omega}{2}\int_{-1}^{1} I(r,\mu') d\mu' + Q(r) \tag{1}
$$

with the boundary condition for µ ∈ [0, 1]:

$$
I(R,-\mu) = K(\mu) + \alpha\, I(R,\mu) + \beta\,\chi(\mu)\int_{0}^{1} I(R,\mu')\mu'\,d\mu' \tag{2}
$$

with `∫₀¹ χ(µ)µ dµ = 1` (Eq. 3).

**Notation map to ORPHEUS** (ORPHEUS = `orpheus.derivations.continuous.peierls`):

| PS-1982     | ORPHEUS                                | Notes                                                   |
| ----------- | -------------------------------------- | ------------------------------------------------------- |
| `r`         | `r` (already in optical units!)        | `r_phys * Σ_t` — code stores cm; PS uses MFP            |
| `R`         | `R · Σ_t`                              | Sphere radius in MFP                                    |
| `µ`         | `µ`                                    | Same convention: cosine wrt outward radial direction    |
| `ω`         | `c = (Σ_s + νΣ_f)/Σ_t`                 | Single-scattering albedo. PS uses 1G transfer-equation  |
| `Q(r)`      | `S(r)/Σ_t`                             | Isotropic internal source per unit optical length       |
| `α`         | (specular BC strength)                 | Variant α takes α=0 (vacuum) or α=1 (white/specular)    |
| `β`         | (diffuse BC strength)                  | Variant α currently has no β                            |
| `K(µ)`      | (incident illumination)                | Variant α currently has no incident source              |
| `I(r,µ)`    | `ψ(r,µ)`                               | Specific intensity = angular flux                       |
| `I(r)`      | `φ(r)`                                 | Energy density = scalar flux (Eq. 7)                    |

### The kernel (Eq. 21) — the central deliverable

PS-1982's final integral equation for `rI(r) = rφ(r)`:

$$
r I(r) = r G(r) + \int_{0}^{R} x\!\left[\tfrac{\omega}{2} I(x) + Q(x)\right]
        \!\Big[E_1(|r-x|) - E_1(r+x) + \alpha F_1(r,x) + \beta F_2(r,x)\Big] dx \tag{21}
$$

The "vacuum source-driven" reduction (α=0, β=0, K=0) collapses to:

$$
r\phi(r) = \int_{0}^{R} x\!\left[\tfrac{\omega}{2}\phi(x) + Q(x)\right]
            \big[E_1(|r-x|) - E_1(r+x)\big] dx \quad \text{(vacuum BC)}
$$

**This is the L1 reference for Plan A2.** Identifying the Peierls
kernel with the standard form: the symmetric-difference part
`[E_1(|r−x|) − E_1(r+x)]` is exactly what falls out of any sphere
Peierls reduction. The transition from `rI(r)` ↔ `I(r)` is the
spherical "rφ" device — multiplying by r symmetrises the radial
problem to a slab on r ∈ [−R, R] (this is the pseudo-slab equivalence
PS cite as the technique of Mitsis 1963 / Erdmann-Siewert 1968 /
Wu-Siewert 1975).

### BC-correction kernels (Eqs. 23–25, vacuum limit α=β=0 → vanish)

The α-branch needed for Plan A2 specular extension:

$$
F_1(r,x) = \begin{cases}
4x \int_{0}^{1} T[\mu_0(x,\mu)]\, \pi^{-1}(r,x,\mu)\cosh(x\mu)\cosh[\pi(r,x,\mu)]
           e^{-2R\mu_0(x,\mu)} d\mu, & x \le r \\[4pt]
4r \int_{0}^{1} T[\mu_0(r,\mu)]\, \pi^{-1}(x,r,\mu)\cosh(r\mu)\cosh[\pi(x,r,\mu)]
           e^{-2R\mu_0(r,\mu)} d\mu, & x \ge r
\end{cases} \tag{23,24}
$$

with the **specular reflection amplifier**:
$$
T(\mu) = [1 - \alpha e^{-2R\mu}]^{-1} \tag{14}
$$
which is the geometric series `Σ α^n e^{-2nRµ}` (n bounces, each
bounce contributing one round-trip optical depth `2Rµ`).
**T(µ) is the integral-operator analog of the
`(1 − α e^{−2τµ})^{−1}` continuum-µ specular kernel that appears in
Sanchez 1986 Eq. (A6) and in Variant α Phase 5.** Sanity check on
the Variant α prototype: `T(µ)` row-summed at α=1 should reproduce
the white-BC closure as a Neumann series resummation.

### Auxiliary kernel functions

- `µ_0(r,µ) = cos θ_0 = [1 − (r/R)²(1 − µ²)]^{1/2}` (Eq. 5) — direction
  cosine on the ray as it hits the surface.
- `S_0(r,µ) = rµ + [R² − r²(1 − µ²)]^{1/2}` (Eq. 6) — chord length
  from (r,µ) to the surface.
- `π(x,r,µ) = [x² − r²(1 − µ²)]^{1/2}` (Eq. 10) — radial component on
  inner integral; vanishes when x = r·sin(θ).

### Mathematical path — structural-independence judgment

**PS-1982 derives the kernel by the "radial-µ + integration" path:**
1. Eq. (4): integrate the integro-differential equation along the
   propagation ray, using the Case-de-Hoffmann-Placzek geometry.
2. Eqs. (8a)+(8b): split into outgoing (µ ≥ 0) and incoming (µ ≤ 0)
   half-spaces, expressed as line integrals with limits
   `r√(1−µ²) → r → R`.
3. Eq. (11): integrate Eqs. (8) over µ and **add the two halves**.
   The `[E_1(|r−x|) − E_1(r+x)]` structure emerges from the
   integration and the explicit chord limits — NOT from any
   `cosh(ρµ)` even-extension argument.
4. Eqs. (12)–(20): solve for the unknown surface flux
   `I(R, −µ)` self-consistently using Marshak-projected balance,
   producing the L(x) kernel (Eq. 20) and J (Eq. 19).
5. Eq. (21): substitute back to close the integral equation.

**Sanchez 1986 (NSE 92) by contrast:** uses the cosh-even-extension
argument that maps a sphere problem to a slab on [−R,R] with a
"reflected" coordinate. Eq. (5) of Sanchez 1986 gives
`ḡ_2(ρ' → ρ) = (ρ'/ρ)·[E_1(σ|r−r'|) − E_1(σ(r+r'))]/2` directly
from cosh-symmetrisation.

**Verdict: STRUCTURALLY INDEPENDENT (verified by direct comparison).**
The two papers arrive at the same `[E_1(|r−x|) − E_1(r+x)]` kernel
by genuinely different mathematical paths — the integration-over-µ
route (PS) vs the slab-extension route (Sanchez 1986). At the
factor-of-2 level the two forms agree (Sanchez has `(ρ'/ρ)·[…]/2`
because it's the angular kernel, PS has `x·[…]` after multiplying
by `rφ(r)`). **This rules out reference contamination of Variant α
α=0 branch.**

PS-1982 explicitly says (last sentence of §2): *"the method of
characteristics has been used to provide an independent verification
of Eq. (21)."* — a third structurally-independent confirmation.

### What PS-1982 does NOT provide

- **No k-eigenvalue values.** Eq. (1) is a fixed-source equation
  with absorption (`ω` enters as scattering ratio). To extract a
  k_eigenvalue Plan A2 must either:
  (a) replace `ω/2` with `c/(2k)` and iterate, or
  (b) interpret the critical condition as the value of `ω` for which
  the homogeneous (Q=K=0) kernel has eigenvalue 1.
- **No numerical method, no test cases, no benchmark values.** PS-1982
  is a derivation-only Note. The numerical work using this kernel is
  what Mitsis 1963 (ANL-6787, in the literature directory) does for
  vacuum sphere — that's the actual numerical-reference paper.
- **No mention of multi-region or multi-group.** Single-region,
  monoenergetic only.
- **No critical-radius tables.** None reported anywhere in the paper.

### V&V framework classification (per vv-principles SKILL.md §"three pillars")

**Pillar:** Semi-analytical reference. Eq. (21) is a Fredholm-2
integral equation; solving it requires numerical quadrature
(Gauss-Legendre on the kernel + Nyström or Galerkin discretisation),
not a closed-form answer. The integrand involves only E_1, cosh, exp,
Gauss-Legendre quadrature on µ ∈ [0,1] — all evaluable to arbitrary
precision via `mpmath`.

**Claim level supported (per the hierarchical claim taxonomy):**
- ✓ **Convergence-order claims:** ψ(r)·r residual converges as
  Plan A2 refines integration mesh.
- ✓ **Flux-shape claims:** the φ(r) profile from Eq. (21) by
  high-order Gauss-Legendre + Nyström is a structurally independent
  reference for Variant α flux-shape verification.
- ✗ **Eigenvalue claims:** NOT directly. PS-1982 is fixed-source.
  However, the **homogeneous-kernel-eigenvalue trick** lifts it:
  finding the value of ω (= c) such that the kernel
  `K(r,x) = (ω/2)·x·[E_1(|r−x|) − E_1(r+x)]` has unit eigenvalue
  on r ∈ (0, R] gives the critical scattering ratio `c_crit(R)`,
  which IS a structurally-independent k-type benchmark for Plan A2.
  This is the trick to use. Plan A2 prototype should perform a
  power iteration on the integral kernel `(ω/2)·x·[E_1 − E_1]`
  with c = ω free, find the largest eigenvalue λ(c), and check that
  λ(c) = 1 yields the critical c for the chosen R. Variant α at
  the same R, c should agree.

### How Plan A2 should use PS-1982

1. Implement the kernel `K(r,x) = x·[E_1(|r−x|) − E_1(r+x)]`
   (vacuum BC, α=β=K=Q=0) using `mpmath.expint(1, ...)` for E_1.
2. Discretise on Gauss-Legendre nodes r_i ∈ (0, R], 50–80 points.
3. Form the matrix `A_ij = w_j · K(r_i, x_j)`.
4. Power iteration on `(c/2)·A` → dominant eigenvalue λ(c).
5. The critical c is `c* = 2/λ_max(A)` where λ_max(A) is the
   largest eigenvalue of A.
6. Cross-check Variant α prototype at (R, c*) — it should give
   k_eff = 1 to many significant figures.
7. Sanity for high mpmath precision (50+ dps): the converged c*
   must agree with classical critical-sphere values from
   Lewis-Miller, Bell-Glasstone, etc. in the c → 1 limit (very
   thick scatterer).

### Suggested L1 reference values to compute (Plan A2 deliverables)

| R (MFP) | c (scattering ratio) | Expected status         | Use as           |
| ------- | -------------------- | ----------------------- | ---------------- |
| 1.0     | various c < 1        | Subcritical, fixed-Q    | Flux-shape L1    |
| various | c = c_crit(R)        | Critical (k_eff = 1)    | Eigenvalue L1    |
| 5.0     | c → 1                | Asymptotic Milne limit  | Asymptotic check |

The c_crit(R) values for the bare critical sphere are tabulated in
Sood-Forster-Parsons 2003 (LANL critical-benchmark report) — that's
the **independent eigenvalue ground** Plan A2 should cite alongside
PS-1982. Plan 2 verification chain is then:
**Variant α (integral op + bounce sum, α=0) ↔ PS-1982 Eq. (21) Nyström
↔ Sood 2003 c_crit(R) tabulated values.**

---

## Part 2 — Garcia 2021 (priority for Plan B3, with caveat)

### Citation
Garcia, R.D.M., "Accurate spherical harmonics solutions for neutron
transport problems in multi-region spherical geometry," J. Comput.
Phys. **424**, 109856 (2021). DOI `10.1016/j.jcp.2020.109856`.

### Problem statement (Eq. 1, abridged)

One-speed transport in a K-region sphere with anisotropic scattering:

$$
\mu \frac{\partial}{\partial r}\Psi_k + \frac{1-\mu^2}{r}\frac{\partial}{\partial\mu}\Psi_k
   + \sigma_{t,k}\Psi_k
   = \frac{\sigma_{s,k}}{2}\sum_{l=0}^{L_k}(2l+1) f_{k,l} P_l(\mu) \int_{-1}^{1} P_l(\mu')\Psi_k d\mu'
   + \frac{\nu_k \sigma_{f,k}}{2}\int_{-1}^{1}\Psi_k d\mu' + S_k(r,\mu)
$$

for k = 1, 2, …, K. **Subcritical: c_k = (σ_{s,k} + ν_k σ_{f,k})/σ_{t,k} < 1.**

### CRITICAL caveat for Plan B3 — no k-eigenvalue results published

**Garcia 2021 page 3, after Eq. (7), states verbatim:**
> *"in this work we focus our attention on the subcritical case
>  c_k < 1 for k = 1, 2, ..., K. Criticality problems will be
>  addressed in a future work."*

**Implication for Plan B3:** This paper does NOT provide reference
k_eff values for the multi-region sphere benchmark. What it does
provide is **converged scalar flux φ_k(r) and current J_k(r) tables
for fixed-source problems** at very high P_N order (up to N=199).
Plan B3 needs to be reframed: the cross-check Garcia 2021 enables
is **scalar-flux and current shape**, not eigenvalue. If a true
multi-region eigenvalue cross-check is wanted, the dedicated
follow-up paper (when it appears, or earlier ones like Williams 2005
Ann. Nucl. Energy 32, Garcia Ref [16] which the paper says was
"interesting" for varying-order P_N but unexplored) must be located
separately.

### Benchmark configuration (the only one in the paper)

3-region sphere from Williams 1991 (Ann. Nucl. Energy 18, 371),
parameters from Garcia 2021 Table 1:

| Region k | R_k (cm) | σ_s,k (cm⁻¹) | σ_t,k (cm⁻¹) | c_k    | a_k (mfp) | b_k (mfp) |
| -------- | -------- | ------------ | ------------ | ------ | --------- | --------- |
| 1 (core) | 3.0      | 0.99         | 1.0          | 0.99   | 0.0       | 3.0       |
| 2 (mid)  | 5.0      | 0.3          | 0.5          | 0.6    | 1.5       | 2.5       |
| 3 (out)  | 7.0      | 1.9          | 2.0          | 0.95   | 10.0      | 14.0      |

Three sub-cases (Garcia 2021 Table 2):
- **Case 1**: internal sources S_1 = 0.5, S_2 = 1.0, S_3 = 1.5;
  no incident flux. Originally Williams 1991 Example 5.
- **Case 2**: no internal sources, incident isotropic flux f(µ) = 1.0
  on the outer surface. Strong attenuation expected through region 3.
- **Case 3**: same as Case 1 but with anisotropic scattering in
  regions 1 and 3 (Legendre coefficients of the scattering law from
  Oblow et al. 1974, NSE 54 — elastic 14.1-MeV n on ²³⁸U for region 1,
  elastic 14.5-MeV n on ²⁰⁸Pb for region 3, with L_k = 14).

### Reference values usable as L1 evidence (NO k_eff)

**Case 1 converged scalar flux** (Garcia 2021 Table 5, rightmost column):

| r (cm) | Region | Converged φ (ppP_N) |
| ------ | ------ | ------------------- |
| 0.0    | 1      | 18.860              |
| 0.5    | 1      | 18.756              |
| 1.0    | 1      | 18.442              |
| 1.5    | 1      | 17.911              |
| 2.0    | 1      | 17.145              |
| 2.5    | 1      | 16.095              |
| 3.0    | 1/2 i. | 14.381              |
| 3.5    | 2      | 13.455              |
| 4.0    | 2      | 13.337              |
| 4.5    | 2      | 13.590              |
| 5.0    | 2/3 i. | 14.361              |
| 5.5    | 3      | 15.532              |
| 6.0    | 3      | 14.198              |
| 6.5    | 3      | 10.807              |
| 7.0    | 3 surf | 4.0763              |

**Case 1 converged current** (Garcia 2021 Table 6 rightmost column):
zero at r=0 (sphere centre, by symmetry), changes sign at r≈4.0 cm
(net inflow to core from boundary illumination not present here →
this is the natural flux-shape pattern).

**Case 2 converged scalar flux** (Garcia 2021 Table 12 rightmost):
external incident flux f(µ) = 1.0; flux maximum at outer surface
(1.6910 at r=7.0 cm), strongly attenuated to 0.26164 at the centre.

**Case 3 converged scalar flux** (Garcia 2021 Table 18 rightmost):
isotropic vs anisotropic comparison; centre flux drops from 18.860
(Case 1) to 12.463 (Case 3) due to forward-peaked scattering
increasing leakage.

**Convergence of P_N method** (Garcia 2021 Table 7):
Case 1 percent RMS errors of P_N flux vs converged ppP_N are 5.3, 1.9,
0.65, 0.32, 0.067 at N = 1, 3, 9, 19, 59. Convergence rate is
sub-exponential (geometric) — a high-order method but not spectral
in N. For Plan B3 cross-checks, one *probably* wants ≤ 5e-4 relative
agreement to call it a match — that means using the rightmost column
("Converged" ppP_N) and not the N=199 P_N column.

### Williams 1991 cross-validated independently

Garcia 2021 Table 8 cross-checks Case 1 against:
- Williams 1991 (integral-equation method via method-of-characteristics
  averaged over shells; same case as proposed it).
- Picca-Furfaro-Ganapol 2012 (NSE 170): discrete-ordinates with
  multiproblem strategy and convergence acceleration.
- Garcia's own ppP_N method.

At most r-positions, three methods agree to 3 significant figures
(reported max difference: 4-th significant figure agreement
everywhere except r=5.0 cm where Williams' and Picca's spatial
discretisation was apparently not converged — Williams later
re-ran with refined grids and confirmed Garcia's 4.0763 at r=7).

**This is a real structural-independence chain.** Three methods
(integral equation MOC + discrete ordinates + spectral P_N), all
agreeing on Case 1 scalar flux. Plan B3 inherits that triangle —
and any future multi-region Variant α extension at this 3-region
configuration becomes cross-checkable against the converged ppP_N
column in Tables 5, 12, 18.

### Stable-P_N machinery (what "stable" means)

The standard P_N solution for spheres has been known since Davison
(1958). The numerical instability comes from condition-number blowup
when the Bessel coefficients are computed in the standard
(non-rescaled) form. Garcia's stabilisation has two levers:

1. **For the central region (Garcia-Siewert-Thomas 2017, Ref [14])**:
   express the modified spherical Bessel functions of the first kind
   using a **scaled ratio** `ι_n(x:y) = i_n(x)/i_0(y)` with `x ≤ y`,
   computed by **backward recurrence**. The ratios stay between 0
   and 1 even when the individual `i_n` overflow. For final
   coefficient determination, use **least-squares + SVD** on the
   linear system from Marshak boundary projection, since the system
   is itself ill-conditioned. This gives the central-region solution
   in ppP_N form.

2. **For the shell regions (Garcia 2019, JCP 405, Ref [15])**:
   reformulate the angular redistribution `(1-µ²)/r · ∂Ψ/∂µ` as a
   **heterogeneous source** in the homogeneous P_N problem, so the
   streaming part becomes plane-geometry-like. The particular
   solution is then obtained by solving an `N+1`-dimensional system
   of **Volterra integral equations of the second kind** by
   Gaussian elimination on the discretised system. This avoids
   the BesselK third-kind ill-conditioning that defeats the
   standard shell P_N solution.

3. **Iterative region-coupling** (this paper, Garcia 2021): sweep
   forward and backward across regions, updating each region's
   incident-flux estimate from the neighbour's most-recent solution.
   Convergence: 10-40 sweeps to reach 10⁻⁹ relative current change.

4. **Post-processing (ppP_N)**: instead of evaluating the angular
   flux from its truncated Legendre expansion (which gives Gibbs-like
   artifacts at angular discontinuities), evaluate via integration
   along the neutron path (Eq. 33, 34) using the previously-computed
   Legendre **moments** to drive the source. This gives the angular
   discontinuities at `µ_c,k(ρ) = [1 − (a_k/ρ)²]^{1/2}` (the
   shadowing critical cosine, Eq. 32) accurately. The current is
   computed from a **flux-conservation balance** on a spherical shell
   (Eq. 40) instead of from the angular-flux integral, since the
   latter loses precision when scattering dominates and the partial
   currents are nearly equal-and-opposite.

### Boundary conditions covered

Garcia 2021 covers: **vacuum (no-incoming)**, **isotropic incident
illumination** `f(µ) = const` (Case 2), and continuity at internal
interfaces. **No specular BC**, **no white BC**, **no reflective
BC** in any reported test case. The method admits incident `f(µ)`
of arbitrary functional form — so reflective/white BC could in
principle be implemented as a self-consistent fixed-point iteration
for `f(µ)` matching the outgoing-current shape — but the paper does
not do this and does not report a numerical example.

### V&V framework classification

**Pillar:** Semi-analytical reference. The ppP_N solution is exact
to spectral (P_N) precision in the angular variable when N is large,
and uses arbitrary-precision arithmetic for the integration along
the neutron path (Volterra solve + post-processing integral). The
correctness ladder (vv-principles §"Semi-analytical correctness ladder"):
1. Integrator correctness: standard SVD/LU/Gauss-Legendre — assumed
   correct (well-tested upstream linear-algebra and quadrature).
2. Reduction correctness: angular-redistribution-as-source recasting,
   scaled-Bessel ratios, Marshak projection — these have been
   independently verified against Williams 1991 and Picca 2012 in
   Garcia 2021 Table 8 to 4 significant figures.

**Claim level supported:**
- ✓ **Convergence-order:** monitor ppP_N residual as N → ∞.
- ✓ **Flux-shape:** Tables 5, 12, 18 are 4-figure references for
  three multi-region scalar-flux profiles.
- ✗ **Eigenvalue:** explicitly *out of scope* for this paper.

**Verdict for Plan B3:** Garcia 2021 is the right reference for
**multi-region scalar-flux shape verification** in a 3-region
fixed-source (Case 1) or external-illumination (Case 2) configuration.
It is NOT the right reference for k-eigenvalue cross-check. Plan B3
should be **rescoped** to use Garcia 2021 for flux-shape evidence
and to source a **separate** multi-region critical-sphere benchmark
(Sood 2003 covers some, Williams 2005 Ann. Nucl. Energy 32 covers
the annular-gap geometry for varying P_N order) for any eigenvalue
claim.

---

## Part 3 — Garcia 2017 / 2019 / 2019b (context)

### Garcia 2017 (JQSRT 196, 155-158) — PARTICULAR SOLUTION
4-pp Note. Builds the *particular* P_N solution `U_j(r), V_j(r)` for
a position-and-angle-dependent source `S(r,µ)`. Key result: Eqs. (13a,b)
+ condition (15) reduce the problem to determining `S_α,n(r)` (the
source moments) and integrating an N+1 system of ODEs in r. Cited in
Garcia 2021 as Ref [21] — the building block when Q_k is not constant
in a region (e.g. when post-processing introduces a position-dependent
collided source).

### Garcia 2019 (JCP 405, 109139) — SPHERICAL SHELL stable solution
Establishes the stable P_N solution for a single annular region
(spherical shell with inner radius a, outer radius b). Solves the
known instability of the standard shell P_N (which uses both first-
and third-kind modified spherical Bessel functions). Reformulates
the angular redistribution as a source so the homogeneous part is
plane-geometry-like; the particular solution requires solving N+1
coupled Volterra integral equations of the second kind by Gaussian
elimination. Reports numerical results for reflection and transmission
probabilities of a shell. Cited in Garcia 2021 as Ref [15] — the
shell solution reused for all intermediate and outer regions of the
multi-region method.

### Garcia 2019b (JCTT 47, 400-423) — EXTERIOR OF SPHERE stable solution
Demonstrates that ill-conditioning in shell P_N comes from the
modified spherical Bessel functions of the **third kind** (k_n),
not the first kind (i_n). Solves the simpler problem (exterior of a
sphere — single shell extending to infinity) where the solution is
expressible in third-kind functions only, and develops a stable
plane-geometry-like recasting. This exterior-sphere result is the
methodological prototype for the shell solution of Garcia 2019, which
is the building block of Garcia 2021. Cited in Garcia 2021 as Ref [38].

### Citation chain for Plan B3 / multi-region work

```
Garcia 2017 (particular P_N)
                    ↘
Garcia-Siewert-Thomas 2017 (stable sphere P_N)
                    ↘
Garcia 2019b (stable exterior P_N → identifies k_n as culprit)
                    ↘
Garcia 2019 (stable shell P_N via Volterra recasting)
                    ↘
Garcia 2021 (multi-region sphere by sweeping)
```

For Plan B3 the load-bearing references are **Garcia 2021** (the
benchmark provider) + **Garcia 2019** (the underlying shell-solution
mechanism, useful if a reproducer of any single sub-case is needed).
Garcia 2017 is only relevant if Plan B3 needs to handle a
position-dependent source (i.e. a previously-collided source from a
preceding sweep).

---

## Cross-paper structural-independence summary

| Reference          | Method                                             | Independent of Variant α? | Independent of Sanchez 1986? | Provides k_eff? |
| ------------------ | -------------------------------------------------- | ------------------------- | ---------------------------- | --------------- |
| PS-1982            | Integral eq. via radial-µ integration              | ✓ STRUCTURAL             | ✓ STRUCTURAL                | ✗ (homog-eigval trick) |
| Sanchez 1986       | Integral eq. via cosh-even-extension               | ✓ STRUCTURAL             | (self)                       | ✗ (fixed-source) |
| Mitsis 1963 (ANL-6787) | Integral eq. + numerical (pseudo-slab)         | ✓ STRUCTURAL             | ✓ STRUCTURAL                | ✓ (critical sphere) |
| Garcia 2021        | Stable spectral P_N + region-sweep                 | ✓ STRUCTURAL             | ✓ STRUCTURAL                | ✗ (subcritical only) |
| Sood 2003 (LANL)   | Compiled critical-sphere benchmarks               | ✓                         | ✓                            | ✓ (tabulated)   |

The **L1 eigenvalue ground for Plan A2 vacuum sphere** is the
combination **PS-1982 Eq. (21) homogeneous-eigenvalue trick + Sood 2003
tabulated c_crit(R)**. Variant α prototype with α=0 should reproduce
both.

The **L1 flux-shape ground for Plan B3 multi-region** is **Garcia 2021
Tables 5, 12, 18 converged ppP_N column**. This is fixed-source, NOT
eigenvalue. Plan B3 must be rescoped accordingly.

---

## Recommended MEMORY.md entry

```
- [PS-1982 + Garcia 2021 extraction](ps1982_and_garcia_extraction.md) — PS-1982 Eq.(21) is structurally-indep vacuum-sphere kernel for Plan A2 (NO k_eff in paper; use homog-eigenvalue trick + Sood 2003 c_crit). Garcia 2021 is multi-region fixed-source ONLY (subcritical c_k<1, criticality "future work") — provides flux-shape Tables 5/12/18 for Plan B3, NOT eigenvalues. T(µ)=[1−α exp(−2Rµ)]^{−1} in PS Eq.(14) is the integral-op analog of Sanchez 1986 Eq.(A6) Phase 5 specular kernel.
```

---

## Action items (literature-researcher recommendations to user)

1. **Plan B3 needs rescoping** before implementation: Garcia 2021 does
   not give k_eff. Either find Garcia's "future work" follow-up paper
   (search OpenAlex / Semantic Scholar for *Garcia ann nucl energy
   multi-region critical sphere* with year ≥ 2021) or use a different
   multi-region eigenvalue benchmark (Williams 2005 ANE 32; ICSBEP-IRPhE
   compilations).

2. **PS-1982 homog-eigenvalue trick is the right Plan A2 implementation**
   — power iteration on the kernel `(c/2)·x·[E_1(|r−x|) − E_1(r+x)]`
   gives c_crit(R) at arbitrary mpmath precision. Cross-check against
   Sood 2003 LA-13511 tabulated values.

3. **The α-branch of PS-1982 Eq. (14)** `T(µ) = [1 − α exp(−2Rµ)]^{−1}`
   is the closed-form integral-operator analog of the continuous-µ
   specular kernel that Phase 5 needs. Sanity check: at α=1 (white BC)
   the Neumann series of T(µ) over µ ∈ [0,1] should reproduce the
   white-BC closure that Phase 5 has been verifying. If Phase 5 gets
   a different answer, the Phase 5 implementation is wrong (or the
   Phase 5 white-BC reference is wrong — but in either case PS-1982
   is the structural-independent oracle).

4. **Suggest adding to Zotero**: Sood-Forster-Parsons "Analytical
   Benchmark Test Set for Criticality Code Verification" (LA-13511,
   2003) — the source of compiled critical-sphere c_crit(R) tables
   that grounds Plan A2 eigenvalue layer. The user should add this
   themselves; the librarian agent does not mutate the user's library.
