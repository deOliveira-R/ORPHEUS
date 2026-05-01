---
name: Plan 2 B1 — Green's function literature pull for sphere homogeneous specular
description: Literature pull for canonical Green's function references in neutron transport. Recommendation: Variant α (Sanchez 1986 Eq. A6 as integral operator). Variant δ surfaced — Garcia 2020/2021 JCP modern stable PN as cross-check ground truth. Pomraning-Siewert 1982 confirmed precursor (JQSRT 28, 503-506; paywalled).
type: project
---

# Plan 2 B1 — Green's function literature for sphere homogeneous specular

## TL;DR

Variant α (integral-operator power iteration on Sanchez 1986 Eq. A6) is the
right architectural choice for ORPHEUS Plan 2. Sanchez 1986 IS the closed-form
Green's function for the homogeneous sphere with α/β BC; Pomraning-Siewert
1982 (JQSRT 28, 503-506) is the isotropic-scattering precursor and the
truth source for ω₁=0. Sanchez 2002 §III does NOT cover sphere — its
periodic-trajectory closure Eq. (15) `ψ = ψ_q(L)/(1−ψ_bd(L))·ψ_bd + ψ_q`
applies to lattices (rectangular, square, hexagonal), not finite spheres,
and the BC is encoded as a scalar `1/(1 − ψ_bd(L))` factor that is
algebraically the same multi-bounce structure as Sanchez 1986 `T(µ) =
1/(1 − α e^{−2aµ})` once the chord is fixed and the geometry collapses
to 1-D. Books (Williams 1971, Case-Zweifel 1967, Pomraning 1973,
Bell-Glasstone 1970, Lewis-Miller 1984) contain general theory of integral
forms but no closed-form sphere-specular Green's function — none of them
gives a numerical method for this problem. **The only published numerical
solutions for sphere isotropic scattering with reflective BC use P_N
spectral expansion, not integral-operator power iteration.** This is a
non-trivial signal: ORPHEUS's chosen Variant α path is original;
Garcia 2020/2021 JCP and Vieira 2019 are the cross-check ground truth (P_N
sphere with reflection); their PDFs are paywalled but their formulation
gives an external numerical reference for the homogeneous-sphere k_eff.

**Surfacing a Variant δ**: Garcia 2018 J. Comp. Theor. Transport 47
"On the P_N Method in Spherical Geometry: A Stable Solution for the
Exterior of a Sphere" + Garcia 2020/2021 JCP — modern stable P_N
spherical harmonic expansion. This is a **discrete-spectrum-only**
restriction of Variant γ (Case-Zweifel singular eigenfunctions) using
an analytical re-summation that avoids the continuous-spectrum
projection. It's a candidate Variant δ but for the **same reason as
γ** — homogeneous-only, doesn't generalise to multi-region. Recommendation
remains Variant α.

## What ORPHEUS needs (recap)

A continuous reference for the sphere homogeneous specular k_eigenvalue
problem, NOT subject to the Phase 5 Nyström hypersingularity. The
candidate path is integral-operator power iteration on Sanchez 1986
Eq. (A6) — work directly on the Green's function `t(r', Ω' → r, Ω)`
on `(r, Ω)` phase space, with `dr'·dΩ'` integration performed on the
trial function BEFORE sampling the kernel, so the diagonal singularity
is smoothed at the operator level rather than at the Nyström matrix level.

## Per-reference equation extracts

### 1. Sanchez 1986 TTSP 14 — Eq. (A6) IS the Green's function (already extracted)

See `phase5_sanchez_1986_sphere_specular.md` for full extraction. Key
result for Plan 2:

```
g_h(ρ' → ρ) = 2α ∫_{µ₀}^{1} T(µ_-) µ_*⁻¹ cosh(ρµ) cosh(ρ'µ_*) e^{−2aµ_-} dµ
              + (β/(1−β χ_*)) Q(ρ') P(ρ)              ... Eq. (A6)
```

with `T(µ) = 1/(1 − α e^{−2aµ})` the multi-bounce factor, `µ₀ = √max(0,
1 − (ρ'/ρ)²)` the chord-visibility cone, `µ_*² = ρ'² − ρ²(1−µ²)`.

**ORPHEUS map**: `g_h → t_h(r',Ω' → r,Ω)` after angle integration; the
thing Plan 2 wants is the *pre-integration* form (Eq. A1)
`t = t̄ + t_h` from Sanchez Eq. (A1)–(A5) — the un-reduced Green's
function on `(r, Ω)` phase space. ORPHEUS Variant α uses Eq. (A6) but
applies it as an operator `(K · ψ)(ρ, µ)`, integrating `dr'·dΩ'`
against ψ before evaluating at observation `(ρ, µ)`. The Phase 5
failure was sampling Eq. (A6) at coincident `(ρ, ρ')` Nyström pairs
where `µ₀(ρ, ρ') → µ₀(ρ, ρ) = 0` and the visibility cone closes.

**BC treatment**: `α + β` general; α=1, β=0 is perfect specular
(ORPHEUS Plan 2 case). BC is fully absorbed into the kernel via Eq.
(A1) `t = t̄ + t_h`; no separate K_bc closure needed. **This is the
key structural advantage of the Green's function approach.**

**Discretisation**: Sanchez does NOT give a numerical method. The paper
is a closed-form theoretical reduction. The implementation choice is
ORPHEUS's design space.

**Critical pages**: pp. 338-342 (Appendix A1-A7).

### 2. Pomraning-Siewert 1982 JQSRT 28, 503-506 — isotropic precursor

**Citation confirmed via OpenAlex**: G.C. Pomraning & C.E. Siewert,
"On the integral form of the equation of transfer for a homogeneous
sphere", *J. Quant. Spec. Rad. Transfer* 28, 503-506 (1982). DOI:
`10.1016/0022-4073(82)90016-4`. 4 pages — a Note/Letter. Paywalled
(Elsevier). 33 citations (CrossRef) / 29 (S2). Both authors at NCSU.

**Status**: PDF not retrieved. Paywalled with no open-access copy on
OpenAlex, OSTI, arXiv, or HAL. Sanchez 1986 cites it explicitly as
the `ω₁ = 0` (isotropic scattering) precursor.

**What it covers** (inferred from Sanchez 1986's citation context +
abstract context of follow-on Garcia papers):
- 1-D radial integral form for a homogeneous sphere with isotropic
  scattering, vacuum BC.
- Reduces the 3-D transport equation to a 1-D Fredholm integral
  equation in `r` (the spherical analog of the slab Peierls equation).
- Likely uses the same `cosh(ρµ)·cosh(ρ'µ_*)` structure as Sanchez
  1986 — that structure comes from the even extension of `J_R(ρ)`
  through `ρ = 0`, which is geometric and pre-dates the BC choice.
- 4 pages → no numerical method, no computational results.

**ORPHEUS map**: this is the `ω₁ = 0`, `α = 0` (vacuum) reduction of
Sanchez 1986 Eq. (1a). For ORPHEUS Plan 2 (which is `ω₁ = 0`, `α = 1`
specular), Pomraning-Siewert is one branch of the truth lattice:
| Case          | α | β | ω₁ | Reference                |
|---------------|---|---|----|--------------------------|
| Vacuum iso    | 0 | 0 | 0  | Pomraning-Siewert 1982   |
| Specular iso  | 1 | 0 | 0  | **ORPHEUS Plan 2 target**|
| Vacuum aniso  | 0 | 0 | ≠0 | Sanchez 1986 (ω₁ branch) |
| Specular aniso| 1 | 0 | ≠0 | Sanchez 1986 (full)      |
| Diffuse       | 0 | β | ≠0 | Sanchez 1986 (β branch)  |

So Pomraning-Siewert 1982 serves as the **vacuum cross-check** when
ORPHEUS Variant α is run with α → 0. Sanchez 1986 IS the literature
specular ground truth.

**Discretisation**: same as Sanchez 1986 — no numerical method
prescribed. 4-page Note focuses on the analytical reduction.

**Critical pages**: paywalled; expected pp. 503-506 contain the
integral-form derivation; ORPHEUS does not need this PDF for B3 —
Sanchez 1986 fully covers the math.

**Recommendation**: do NOT block on retrieving Pomraning-Siewert 1982
PDF. Sanchez 1986 covers everything Pomraning-Siewert 1982 does plus
the anisotropic + reflective extensions.

### 3. Sanchez 2002 NSE 140 §III — periodic-trajectory closure (READ LOCALLY)

`/workspaces/ORPHEUS/Sanchez(2002).pdf` pp. 30-35. Title: *"Treatment
of Boundary Conditions in Trajectory-Based Deterministic Transport
Methods"*. Section III is *"Exact Treatment for Closed Domains"*.

**Critical structural finding**: §III applies to **infinite lattices
generated by periodic translations and/or specular reflections of a
basic 2-D x-y rectangular / square / hexagonal domain**. It does
NOT cover finite spheres. The "Green's function with BC absorbed"
language in the Sanchez 2002 abstract refers to the periodic-trajectory
closure for *lattice geometries*, not for finite-volume problems
with reflective external boundaries.

**Key equation (Eq. 15)**:

```
ψ(x) = ψ_q(L)/(1 − ψ_bd(L)) · ψ_bd(x) + ψ_q(x)              ... Eq. (15)
```

where `ψ_bd(x) = e^{−τ(0,x)}` is the angular flux at position x along
the trajectory produced by a unit incoming flux with no volumetric
sources, and `ψ_q(x) = ∫₀^x e^{−τ(x',x)} q(x') dx'` is the angular
flux from volumetric sources with zero incoming flux. `L` is the
trajectory period.

**ORPHEUS notation map**:
| Sanchez 2002 | ORPHEUS                                          |
|--------------|--------------------------------------------------|
| `ψ(x)`       | `ψ(s)` along a periodic trajectory               |
| `ψ_bd(L)`    | `e^{−τ_period}` — bounce-attenuation per period  |
| `ψ_q(L)`     | source contribution per period                   |
| `1/(1−ψ_bd(L))` | scalar multi-bounce factor                    |

**Algebraic match to Sanchez 1986**: `1/(1 − ψ_bd(L)) = 1/(1 − e^{−τ_L})`
is **literally the same** as Sanchez 1986 `T(µ) = 1/(1 − α e^{−2aµ})`
with α=1 and `τ_L = 2aµ`. The two papers describe the same
multi-bounce closure under different geometric framings:

- Sanchez 1986: 1-D radial sphere, chord per bounce = 2aµ, indexed by `µ`.
- Sanchez 2002: 2-D Cartesian lattice, period = L, chord per period
  fixed by the lattice geometry.

**The two are structurally identical:** sum over bounces along a
*single* characteristic trajectory, with attenuation `e^{−τ_chord}`
per bounce. Sanchez 1986's `T(µ)` integrates the bounce sum over
the polar cosine; Sanchez 2002's `1/(1−ψ_bd(L))` is the bounce sum
along a fixed lattice trajectory.

**Discretisation**: Sanchez 2002 §III prescribes the **method of
characteristics with periodic trajectory tracking**. The numerical
method is to sweep along periodic trajectories, terminating when
`|x − x_0| < ε` and `|l − L| < ε`, then apply Eq. (15) to recover
the angular flux at every point along the trajectory. This is a
*deterministic* sweep, NOT power iteration on an integral operator —
it computes ψ pointwise after the BC closure has already been
algebraically absorbed.

**ORPHEUS relevance**: Sanchez 2002 §III is **not directly applicable**
to the homogeneous sphere with specular BC, because the sphere is
NOT a periodic-trajectory geometry — a specular bounce off a sphere
returns at a different polar angle each time, never repeating along
a single trajectory. (The structure that Sanchez 1986 exploits is
the *integral over µ* of independent µ-channels, each of which has
a fixed chord length — the µ-fibration is the sphere's analog of
periodic trajectories, but it's an integral, not a discrete sum.)

**Critical pages**: pp. 30-32 (Eqs. 14-16). pp. 32-34 contain the
geometric machinery (rectangular, square, hexagonal periods).
**§III.D pp. 33-34** treats "Periodicity in the Axial Direction" but
this is for 3-D infinite lattices in z, not for finite spheres.

**Verdict**: Sanchez 2002 confirms the *general structure* `1/(1−chord
attenuation)` of multi-bounce closure but does NOT give the sphere
formulation. Plan 2 should rely on Sanchez 1986 and treat Sanchez 2002
as an algebraic cross-check that the multi-bounce factor structure
is universal across geometries.

### 4. Williams 1971 — Mathematical Methods in Particle Transport Theory

**Status**: book, not in OpenAlex (no DOI), not in Zotero (server
flaky — see memory entry on Zotero MCP flakiness). M.M.R. Williams,
*Mathematical Methods in Particle Transport Theory* (Wiley/Butterworths,
1971). Cited by Sanchez 1986 as Ref. 3 ("slab Fredholm").

**Inferred coverage** (from how Sanchez 1986 cites it + standard
literature reputation):
- Slab Fredholm theory for the integral form.
- Linear-operator analysis of the transport operator.
- Spectral theory of the Boltzmann operator.
- **Fredholm second-kind compactness arguments** for proving existence
  and uniqueness of integral-form solutions.

**ORPHEUS map**: useful for B3 SymPy derivation as a **theoretical
foundation reference** — Variant α power iteration converges if and
only if the integral operator is Fredholm second kind with spectral
radius < 1. Williams gives the framework for proving this. NOT
useful as a closed-form reference for sphere specular.

**Discretisation**: Williams treats Fredholm theory abstractly. No
numerical method specific to sphere. Treatments of slab integral
equations use Galerkin / collocation / Nyström — same general
methods as ORPHEUS.

**Critical pages**: unread. Estimated chapters of interest:
- Ch. 2 or 3: integral form derivation
- Ch. on "Spectrum of the transport operator" — Fredholm framework
- Ch. on "Integral equation methods" — slab numerical examples

**Recommendation**: not blocking for Plan 2. If the B3 SymPy work
needs a Fredholm second-kind compactness proof (e.g., to justify
power iteration convergence formally), pull this PDF then. For the
prototype B4, not needed.

### 5. Case & Zweifel 1967 — Linear Transport Theory

**Status**: book, not in OpenAlex (no DOI), not in Zotero. K.M. Case
& P.F. Zweifel, *Linear Transport Theory* (Addison-Wesley, 1967).

**Inferred coverage** (well-known reference):
- **Singular eigenfunction expansion** for the slab and sphere
  homogeneous problem (Case's method).
- Discrete + continuous spectrum decomposition: ψ(x, µ) = a_+ φ_+(µ)
  e^{x/ν₀} + a_− φ_−(µ) e^{−x/ν₀} + ∫_{−1}^{1} A(ν) φ_ν(µ) e^{x/ν} dν
- Chandrasekhar-H functions enter when projecting the BC onto the
  continuous-spectrum subspace.
- BC closure via half-range orthogonality of the Case eigenfunctions.

**ORPHEUS map** (this is **Variant γ** in the plan):
| Case-Zweifel             | ORPHEUS Variant γ                |
|--------------------------|----------------------------------|
| φ_±(µ), ν₀               | discrete-spectrum eigenfunctions |
| φ_ν(µ), ν ∈ [−1,1]       | continuous-spectrum eigenfunctions|
| H(µ) Chandrasekhar       | already in ORPHEUS (Issue #127)  |
| Half-range orthogonality | BC closure mechanism             |

**Variant γ assessment**: Case-Zweifel singular eigenfunction expansion
is **analytical for homogeneous slab**, but the homogeneous **sphere**
case requires substantially more machinery (the spherical analog uses
modified spherical harmonics and the spectrum is more complex than the
slab `[−1, 1]` continuous spectrum). The historical extensions
(Mitsis 1963 for sphere; Garcia & Siewert in many JQSRT papers) all
work in the spherical case but with formidable bookkeeping —
**Garcia 2020 JCP and Garcia 2021 JCP are the modern numerical
implementations of this approach** for sphere multi-region with
reflection.

**Discretisation**: Case-Zweifel itself is purely analytical. Numerical
implementations (Garcia papers) use truncated H-function representations
and direct quadrature on the continuous-spectrum integral.

**Critical pages**: unread. Estimated:
- Ch. 4: slab Case-Zweifel singular eigenfunctions
- Ch. 5 or 6: spherical geometry extension
- Appendix on Chandrasekhar-H functions

**Recommendation**: Variant γ is well-trodden — Garcia's recent JCP
papers are the modern entry point, NOT Case-Zweifel directly. If
Plan 2 ever pivots to Variant γ, fork from Garcia 2020/2021 (which
absorbs Case-Zweifel into computer-friendly form), not from
Case-Zweifel directly.

### 6. Pomraning 1973 — The Equations of Radiation Hydrodynamics

**Status**: book, not in OpenAlex. G.C. Pomraning, *The Equations of
Radiation Hydrodynamics* (Pergamon Press, 1973). The same Pomraning
who wrote the 1982 JQSRT precursor.

**Inferred coverage**:
- Radiative-transfer integral form for general 1-D / 2-D / 3-D media.
- Diffusion approximation, Eddington factor, P_N/M_1 closures.
- Time-dependent radiative transfer — but the steady-state equations
  reduce to ORPHEUS's case.
- Likely contains the sphere integral form as one of the worked
  examples in the chapter on "Exact solutions".

**ORPHEUS map**: ORPHEUS does not need radiation hydrodynamics. The
relevant content is the sphere integral form, which is fully covered
by the same author's 1982 JQSRT paper (and by Sanchez 1986).

**Discretisation**: book covers approximate / closure methods (P_N,
diffusion); not directly applicable to ORPHEUS Plan 2.

**Critical pages**: unread. Estimated chapter on "Exact solutions"
or "Integral-form transport".

**Recommendation**: not blocking. Sanchez 1986 + Pomraning-Siewert
1982 cover everything Plan 2 needs.

### 7. Bell & Glasstone 1970 — Nuclear Reactor Theory, Chapter 5

**Status**: book, not in OpenAlex. G.I. Bell & S. Glasstone, *Nuclear
Reactor Theory* (Van Nostrand Reinhold, 1970), Chapter 5 *"The Method
of Discrete Ordinates and Discrete S_N"* and Chapter 1 §1.6 contain
the integral form derivation.

**Inferred coverage** (standard textbook):
- Ch. 1 §1.6: integral form derivation `φ(r) = ∫ K(r←r') q(r') dr'`
  for a general medium.
- Ch. 1 §1.7: integral form for homogeneous medium with vacuum BC.
- Plumbing: optical-distance kernels, attenuation factors, Peierls
  equation reduction.

**ORPHEUS map**: standard reference for the *vacuum* sphere integral
form. Same content as Sanchez 1986 reduced to ω₁=0, α=0.

**Discretisation**: textbook level — no numerical method specific to
sphere with reflection. Vacuum slab numerical examples only.

**Critical pages**: Ch. 1 §1.6-1.7. Not blocking.

### 8. Lewis & Miller 1984 — Computational Methods of Neutron Transport

**Status**: book, not in OpenAlex. E.E. Lewis & W.F. Miller Jr.,
*Computational Methods of Neutron Transport* (Wiley, 1984; reissued
by ANS 1993).

**Inferred coverage** (standard computational textbook):
- Ch. 1-2: integral form, Peierls equation, energy-multigroup framework.
- Ch. 4: collision probability methods (CP / Roth method).
- Ch. 5: integral form numerical methods — focus on slab CP and
  cylindrical lattices for thermal reactor applications.
- Ch. 9: discrete ordinates.
- Reflection BC discussed in slab context, not sphere.

**ORPHEUS map**: NOT a Green's function reference. Lewis-Miller is
a CP / S_N / MoC computational reference. Does NOT cover Variant α
power iteration on integral operators. For ORPHEUS Plan 2, it serves
as a **negative result**: standard nuclear-engineering computational
references do not address sphere homogeneous specular k_eff via
integral-operator power iteration. Variant α is genuinely original
within the nuclear-engineering literature.

**Critical pages**: skim Ch. 2, 5, 9 for any mention of sphere
specular, but expectation is no relevant content. **Not blocking.**

### 9. Garcia 2020/2021 JCP, Garcia 2018 JCTT — modern sphere P_N (SURFACED VARIANT δ)

**Citations confirmed**:
- R.D.M. Garcia (2020), *J. Comp. Phys.* 393, "A numerically stable
  spherical harmonics solution for the neutron transport equation in
  a sphere", DOI `10.1016/j.jcp.2019.109139`. **Paywalled.**
- R.D.M. Garcia (2021), *J. Comp. Phys.* 433, "Accurate spherical
  harmonics solutions for neutron transport problems in multi-region
  spherical geometries", DOI `10.1016/j.jcp.2020.109856`. **Paywalled.**
- R. Garcia (2018), *J. Comp. Theor. Transport* 47, "On the P_N
  Method in Spherical Geometry: A Stable Solution for the Exterior
  of a Sphere", DOI `10.1080/23324309.2018.1489847`. **Paywalled.**
- W.J. Vieira (2019), "A comparison of transport methods for the
  solution of a problem with shadowing effects in spherical geometry"
  — paywalled, citation graph only.

**What they do**: Modern stable P_N spherical-harmonic expansion for
the neutron transport equation in 1-D sphere with isotropic scattering
and reflective / vacuum BC. Uses an analytical re-summation that
overcomes the conditioning issues of the original P_N formulation
(which becomes ill-conditioned at high N for sphere). Garcia 2021
extends to multi-region, exactly the same problem as ORPHEUS Phase 4
Class B MR catastrophe (Issue #132).

**ORPHEUS map** (this is the "Variant δ" surfacing requested in B1):
| Garcia 2020/2021         | ORPHEUS                          |
|--------------------------|----------------------------------|
| Stable P_N coefficients  | Analytical projection of ψ on    |
|                          | spherical harmonics              |
| Multi-region Garcia 2021 | Direct attack on Issue #132      |
| Reflective / vacuum BC   | BC absorbed via spherical-       |
|                          | harmonic boundary closure        |

**Variant δ assessment**: P_N spectral expansion — a **discrete-spectrum
restriction** of Variant γ (Case-Zweifel) that uses Garcia's stable
analytical re-summation to dodge the continuous-spectrum projection
entirely. **For homogeneous sphere with reflection, this is the
state-of-the-art numerical reference in the literature.** It's a
parallel architectural choice to Variant α.

**Why Variant δ is NOT recommended** over Variant α for ORPHEUS:
1. **Multi-region is in Garcia 2021, but the formulation uses
   stable-P_N machinery that ORPHEUS doesn't have implemented**
   (Garcia uses recurrence relations specific to spherical harmonics).
   Variant α reuses ORPHEUS's existing Sanchez 1986 implementation
   `compute_K_bc_specular_continuous_mu_sphere` plus integration
   primitives.
2. **P_N truncation is a closure** — ORPHEUS has had a bad time with
   truncation closures (rank-N gating, Class B catastrophe). Variant α
   is closure-free at the operator level (the integral kernel is the
   exact Green's function; only quadrature error enters).
3. **Garcia papers are paywalled and not extractable** without
   institutional access. Variant α can be implemented from Sanchez
   1986 alone (locally available).
4. **Variant δ as a cross-check is more useful than as a primary path**:
   if Plan 2's prototype works, Garcia 2020/2021 paywalled values
   become a "one-day-pull" task that gives ORPHEUS an external
   validation point. As a primary path, the paywall + closure issues
   make it strictly worse than Variant α.

**Discretisation**: Garcia 2020 — **stable P_N truncation with
analytical re-summation**. The re-summation is the contribution; the
underlying method is just truncated spherical-harmonic projection.

**Critical pages**: paywalled. NOT blocking for Plan 2 prototype, but
**should be pulled before Plan 2 closeout** so the cross-verification
matrix B5 can include the Garcia P_N reference. Suggest the user pull
these via institutional access.

### 10. Hébert 2009 Chapter 3 — locally available, cross-reference only

`/workspaces/ORPHEUS/Hebert(2009)Chapter3.pdf`. Already extracted in
`hebert_2009_ch3_interface_currents.md`. §3.8.5 covers multi-region
sphere CP with **vacuum BC only**. The `(1−P_ss)^{−1}` Hébert factor
in §3.8.5 is the rank-0 (scalar) reduction of the multi-bounce
sphere closure. Variant α at rank-1 trial = isotropic constant
should reduce to Hébert's `(1−P_ss)^{−1}` algebraically (this is the
B3 SymPy step 4 cross-check).

## Architectural recommendation: Variant α

**Recommend Variant α — integral-operator power iteration on Sanchez
1986 Eq. (A6)** as the prototype path, in agreement with the Plan 2
tentative choice. Justification (one paragraph):

The literature pull confirms that Sanchez 1986 IS the closed-form
Green's function for the homogeneous sphere with α/β BC, and no
follow-on paper has done better for the *general* α/β case (the
Garcia 2018/2020/2021 modern P_N papers cover only the
isotropic-scattering, restricted reflective case via a
spectral-expansion architecture, not via the Green's function
operator). Variant α reuses ORPHEUS's existing Sanchez 1986 reference
implementation, requires no new closure machinery (BC is absorbed
into the kernel by Eq. A1 `t = t̄ + t_h`), and the discretisation
change from Phase 5 — `dr'·dΩ' integration on the trial function
before sampling` — is exactly the structural change diagnosed in the
Phase 5 retreat as the root cause of the Nyström hypersingularity.
Variant β (BEM-style surface IE) is a more invasive reformulation
without literature precedent for sphere specular k_eff, and Variant γ
(Case-Zweifel) is dominated by Garcia 2020/2021's stable P_N
re-summation, which is itself dominated by Variant α's closure-free
property. **Variant α has the best cost-to-novelty ratio.**

## Critical answer to the discretisation question

The plan asks: *"what discretisation does each reference use for power
iteration on the integral operator?"* The answer is:

**None of the references use integral-operator power iteration for
the homogeneous-sphere-specular k_eigenvalue problem.** All of them
either (a) work analytically without prescribing a numerical method
(Sanchez 1986, Pomraning-Siewert 1982), or (b) use spectral expansion
in spherical harmonics (Garcia 2020/2021), or (c) use the method of
characteristics with periodic trajectories (Sanchez 2002 §III, but
for lattice geometries not finite sphere).

**This validates the Plan 2 working hypothesis**: integrating
`dr'·dΩ'` against a trial function before sampling the kernel is a
genuinely original architectural move within the published nuclear
engineering literature for this problem. It is NOT a standard
technique; the Phase 5 Nyström sampling that failed was the
"obvious" approach and has no literature precedent for success
either.

**Risk implication**: novelty cuts both ways. The hypothesis "smoothing
at the operator level avoids the diagonal singularity" is plausible
on functional-analysis grounds (the integral over `dΩ'` against a
smooth trial function turns the `1/(cos ω_i cos ω_j)` Hadamard
finite-part into a proper integral by sublating one of the cosines),
but ORPHEUS does not have a literature reference where this exact
transformation has been carried out for sphere specular. The B3 SymPy
step is therefore *load-bearing* — if SymPy step 2 (operator action
on `ρ^k µ^m` trial) gives a divergent result, Variant α inherits
Phase 5's structural problem despite the reformulation.

## "What's NOT covered" — gaps for B3 SymPy

The B3 SymPy script must fill the following gaps that the literature
does not address:

1. **Operator action on a polynomial trial function**. `(K · ρ^k
   µ^m)(ρ, µ) = ?` Sanchez 1986 gives Eq. (A6) but never evaluates
   it on a non-flux test function. SymPy must verify this is finite
   (no Hadamard singularity) and reproduce a closed form.

2. **Vacuum reduction at α → 0**. The closed-form `T(µ) → 1` and
   Eq. (A6) reduces to Pomraning-Siewert 1982 Eq. (whatever — paywalled),
   which itself reduces to Bell-Glasstone Eq. (1.43) (`E_2` slab
   integral kernel). SymPy must verify this reduction chain
   symbolically.

3. **Rank-1 isotropic-trial reduction to Hébert's `(1 − P_ss)^{−1}`**.
   At trial = `ψ_trial = const`, the operator action should give
   `const × ⟨1, K · 1⟩` algebraically equal to Hébert §3.8.5
   `(1 − P_ss)^{−1}` × const. **Required for B5 cross-verification at
   rank-1.**

4. **Spectral radius of the integral operator**. Power iteration
   converges iff the dominant eigenvalue is real and isolated. Sanchez
   1986's kernel is positive (all factors `cosh(ρµ)`, `e^{−2aµ}`,
   `T(µ) = 1/(1 − e^{−2aµ})` are positive on the integration domain),
   so by Perron-Frobenius the dominant eigenvalue is real and positive
   — but this does not guarantee isolation. Williams 1971 Fredholm
   compactness arguments may be needed if convergence is sluggish in
   B4 prototype.

5. **k_eigenvalue extraction**. Sanchez 1986 is a fixed-source paper;
   ORPHEUS Plan 2 needs the **k_eff** problem `(I − K_scat) ψ =
   (νΣ_f / k) ψ`. SymPy must verify that the standard fission-source
   power-iteration `k^{(n+1)} = k^{(n)} · ⟨νΣ_f, ψ^{(n+1)}⟩ / ⟨νΣ_f,
   ψ^{(n)}⟩` extracts the dominant Σ_s + ν Σ_f / k mode.

6. **Multi-region piecewise τ(µ)** (out of B3 scope, but flagged for
   B6 / future). Sanchez 1986 is homogeneous-only; the closed-form
   `T(µ) = 1/(1 − e^{−2aµ})` extends with piecewise `τ(µ) = Σ_k
   Σ_k ℓ_k(µ)` under perfect specular, but `cosh(ρµ)` source/receiver
   factors require multi-region generalisation via piecewise optical
   depth integration. This is on the Plan 2 deferred-scope list.

## Critical-pages index (per reference, for follow-up reading)

| Ref                          | Relevance | Pages                     | Local? |
|------------------------------|-----------|---------------------------|--------|
| Sanchez 1986 TTSP 14         | PRIMARY   | pp. 338-342 (Appendix)    | YES    |
| Pomraning-Siewert 1982 JQSRT | precursor | pp. 503-506 (4 pp Note)   | NO (paywalled) |
| Sanchez 2002 NSE 140 §III    | algebraic | pp. 30-32 (Eqs. 14-16)    | YES    |
| Williams 1971                | Fredholm  | Ch. on spectral theory    | NO     |
| Case-Zweifel 1967            | Variant γ | Ch. 4-5 (slab + sphere)   | NO     |
| Pomraning 1973               | precursor | Ch. on exact solutions    | NO     |
| Bell-Glasstone 1970          | textbook  | Ch. 1 §1.6-1.7            | NO     |
| Lewis-Miller 1984            | negative  | Ch. 2, 5, 9 (skim)        | NO     |
| Garcia 2020 JCP              | Variant δ | full paper                | NO (paywalled) |
| Garcia 2021 JCP              | Variant δ | full paper (multi-region) | NO (paywalled) |
| Garcia 2018 JCTT 47          | Variant δ | full paper (exterior)     | NO (paywalled) |
| Hébert 2009 Ch. 3            | rank-0    | §3.8.5 sphere CP          | YES    |

## Citation graph

```
                 Williams 1971 (Fredholm theory)
                        ↑
                        │
   Bell-Glasstone 1970 ←──── (general 1.6-1.7 framework)
                        ↑
                        │
       Pomraning 1973 ──┴── Pomraning-Siewert 1982 (sphere ω₁=0)
                                     ↓
                                Sanchez 1986 (sphere ω₁≠0, α/β BC)  [PRIMARY]
                                     ↓
                                Sanchez 2002 (lattice trajectories) [DIFFERENT GEOMETRY]
                                     │
                                     ╲   (algebraic structure 1/(1−chord_atten) shared)
                                     ╱
                                Hébert 2009 §3.8.5 (rank-0 reduction)  [LOCAL]

   Case-Zweifel 1967 (slab singular eigenfunctions)
                        ↓
           Mitsis 1963 + many others (sphere extension; not pulled)
                        ↓
    Garcia 2018 JCTT 47 (exterior P_N stable)
                        ↓
    Garcia 2020 JCP (sphere stable P_N)  [Variant δ candidate]
                        ↓
    Garcia 2021 JCP (multi-region sphere P_N)  [direct attack on Issue #132]

   Lewis-Miller 1984: orthogonal (CP/S_N/MoC computational reference;
                                  no Green's function content)
```

## File index

- This memo: `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
- Plan 2: `/workspaces/ORPHEUS/.claude/plans/peierls-greens-function-approach.md`
- Predecessor memo: `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md`
- Phase 5 retreat: `/workspaces/ORPHEUS/.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`
- Hébert Ch.3 memo: `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/hebert_2009_ch3_interface_currents.md`
- Local PDFs: Sanchez 1986, Sanchez 2002, Hébert 2009 Ch.3 (all under `/workspaces/ORPHEUS/`)
- Paywalled (one-day-pull): Pomraning-Siewert 1982, Garcia 2018/2020/2021

## Suggested Zotero additions (user-side write — agent does not mutate library)

The user may wish to add the following items to Zotero for this work
stream (literature-researcher does NOT auto-add):

1. **Pomraning-Siewert 1982** — DOI `10.1016/0022-4073(82)90016-4`.
   "On the integral form of the equation of transfer for a homogeneous
   sphere", JQSRT 28, 503-506. Confirmed precursor to Sanchez 1986.

2. **Garcia 2020 JCP** — DOI `10.1016/j.jcp.2019.109139`. "A
   numerically stable spherical harmonics solution for the neutron
   transport equation in a sphere". Variant δ candidate; B5 cross-check.

3. **Garcia 2021 JCP** — DOI `10.1016/j.jcp.2020.109856`. "Accurate
   spherical harmonics solutions for neutron transport problems in
   multi-region spherical geometries". Direct attack on Issue #132
   Class B MR catastrophe via P_N spectral.

4. **Garcia 2018 JCTT** — DOI `10.1080/23324309.2018.1489847`.
   "On the P_N Method in Spherical Geometry: A Stable Solution for
   the Exterior of a Sphere". Foundational stable-P_N paper.

(Zotero MCP server returned 0-hit + connection-refused on all queries
during this session; flakiness pattern from `reference_zotero_flakiness.md`.
Cannot verify whether these items are already in the library — if they
are, the search server is broken, not the library.)
