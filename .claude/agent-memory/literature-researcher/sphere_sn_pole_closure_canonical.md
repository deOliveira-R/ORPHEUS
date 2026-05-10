---
name: Sphere SN pole closure canonical references
description: Hébert §3.9.4 IS the canonical 1D-sphere SN pole stencil; "Bailey 2009" docstring should cite Bailey-Morel-Chang 2010 NSE 165, NOT Bailey-Adams-Yang-Zika 2008 JCP (citation error in ORPHEUS reduced_operator.py)
type: reference
---

# Sphere SN angular-redistribution pole closure — canonical literature

**Cite when**: Issue #168 Defect 3 work, Wave G/H of SN reshape campaign, any
discussion of curvilinear SN apply-vs-sweep operator equivalence at i=0.

## Canonical primary source — LOCAL

**Hébert (2009), *Applied Reactor Physics*, Ch.3 §3.9.4 pp. 141-144**.
File: `scratch/literature/Hebert(2009)Chapter3.pdf`. Contains:

- Eq. (3.418): conservative-form 1D sphere SN, integrated over angular sub-domain
- Eq. (3.420): α-coefficient evaluation of `∫ ∂_μ[(1-μ²)ψ] dμ`
- Eqs. (3.423)-(3.424): α-recursion `α_{1/2}=0`, `α_{n+1/2} = α_{n-1/2} − 2 𝒲_n μ_n`
- Eq. (3.428): cell-balance with redistribution divisor `ΔS_i / (2 𝒲_n)`
  (NOTE the factor of 2 — the apply form `redist = -μ·ΔA[0]·ψ/V[0]` is the
  flat-flux-collapsed form of this, NOT a different equation)
- Eqs. (3.432)-(3.435): **Carlson starting-direction at μ=−1** (the inward zero-weight
  sweep that initializes the α-cascade)
- Eqs. (3.436)-(3.439): asymmetric per-ordinate sweep with DD angular closure

**Lewis & Miller §4.5 NOT NEEDED**: Hébert §3.9.4 contains identical equations.
L&M is paywalled / out-of-print; do not pursue.

## Canonical justification — α-recursion and τ_mm clamp

**Bailey, Morel & Chang (2010), NSE 165(2):149-169** — "Asymptotic Diffusion-Limit
Accuracy of Sn Angular Differencing Schemes," LLNL preprint LLNL-JRNL-420356.
Open-access preprint at https://www.osti.gov/servlets/purl/1020346 .
DOES re-derive the M-M weighted-diamond clamp with formal-ε analysis.
DOES NOT discuss apply-vs-sweep equivalence or i=0 truncation behavior on
non-flat ψ.

## CITATION CORRECTION needed in `orpheus/geometry/reduced_operator.py`

Lines 61-71 cite "Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
*A piecewise linear finite element discretization of the diffusion equation for
arbitrary polyhedral grids*. JCP 227, 3738-3757."

**This is the wrong Bailey paper.** Adams-Yang-Zika is a piecewise-linear FE
DIFFUSION paper — has nothing to do with curvilinear SN α-recursion. The
intended reference is **Bailey-Morel-Chang 2010 NSE 165**. Both papers share
author T. S. Bailey but are unrelated.

## Pomraning 1989 — theoretical justification (LOCAL)

**Pomraning (1989), NSE 101:330-340** "The Transport Equation in General Geometry."
File: `scratch/literature/Pomraning(1989)The Transport Equation in General Geometry.pdf`.
Page 339 explicitly states the sphere center is a structurally-singular point of
the streaming operator: `1/r` factors are intrinsic, "the attendant difficulties
are well known, particularly in numerical treatments." Use this paper as the
"why a separate pole stencil is architecturally justified" reference.

## NOT relevant to Defect 3

- **Morel 1989 NSE 101:72-87** "Hybrid Collocation-Galerkin-Sn" (LOCAL): treats
  scattering source moment-to-discrete matrix only. Page 5 explicitly endorses
  "the standard S_n method should be used for the angular leakage term."
- **Stacey 2007 Ch.9** (LOCAL): Cartesian + integral transport + CP only. Table
  9.1 p.309 has continuum streaming operators (cross-check value); no discrete
  SN sphere.

## Architectural framing

The pole stencil is in the **angular** half of the operator, NOT the spatial
half. Defects 1 and 2 of Issue #168 are spatial (a `BoundaryFaceFlux` Protocol
fits). Defect 3 is angular and needs a separate Protocol — `PoleAngularClosure`
or similar — with two implementations: legacy `BaileyFlatFluxRedist` (the bug)
and canonical `MorelMontryAngularSweep` (Hébert Eqs. 3.432-3.439).

The closure works because the SWEEP processes ordinates μ=−1 → μ=+1 with
asymmetric DD closure `ψ_{n+1/2} = 2ψ_n − ψ_{n-1/2}`; at μ=−1 the redistribution
drops out because `α_{1/2}=0`, giving a slab-like equation that initializes the
α-cascade. The APPLY matvec must reproduce this fixed point; the symmetric
average closure `ψ_{n+1/2} = (ψ_n + ψ_{n+1})/2` does NOT, by exactly the
factor-of-2 the design memo identifies.
