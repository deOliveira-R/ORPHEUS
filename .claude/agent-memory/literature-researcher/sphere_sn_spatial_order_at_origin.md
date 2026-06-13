---
name: Sphere SN spatial convergence order at r=0 (central cell)
description: The literature on O(h^p) at the spherical SN central cell. Bottom line — NO canonical reference reports an O(h^2) central-cell SPATIAL closure; Hébert/Lathrop/Bailey all study ANGULAR error and use uniform DD spatially. First-order-at-origin is an under-documented mesh effect, not a published fundamental result. Cite for any "is O(h^2) at r=0 achievable" question.
type: reference
---

# Spherical SN — spatial convergence order at the r=0 central cell

**Cite when**: any question about pointwise/central-cell SPATIAL convergence
order (O(h^p)) at the sphere coordinate singularity r=0; the ORPHEUS
"O(h^0.95) at innermost cell, O(h^2) interior" diagnosis; choosing a target
reference closure for a second-order-at-origin fix.

## THE decisive disambiguation (notation trap)

Three different "orders" collide in this literature. Keep them apart:

1. **Spatial mesh order O(h^p)** — convergence of the discrete solution to
   the continuous one as Δr→0 at FIXED N. THIS is the ORPHEUS question.
2. **Angular order (in N)** — Lathrop's "quadratic"/"fourth-order error
   reduction" as N→∞ at FIXED/zero spatial error. NOT the ORPHEUS question.
3. **Asymptotic diffusion-limit order (ε-expansion)** — Bailey-Morel-Chang's
   "leading order" vs "first order" preservation of the diffusion limit by
   the ANGULAR scheme. NOT a mesh order. NOT the ORPHEUS question.

The canonical "sphere SN" papers all study (2) or (3), NOT (1). This is why
the literature does NOT contain a published "O(h^2) central-cell spatial
closure" — the papers deliberately remove spatial differencing.

## Hébert (2009) §3.9.4 — LOCAL, the canonical spatial closure

`scratch/literature/Hebert(2009)Chapter3.pdf` pp. 141-144, Eqs.
(3.418)-(3.439). VERIFIED by direct PDF read. Hébert's spatial scheme at
r=0 is just the STANDARD diamond difference (DD), Eq. (3.431)
`φ_{n,i} = ½(φ_{n,i-1/2}+φ_{n,i+1/2})`, applied uniformly. The geometric
weights are V_i, S_{i±1/2} from Eq. (3.427) with `r_{1/2}=0` (explicit:
"with r_{1/2}=0"), so the innermost inner-face area S_{1/2}=4π·0²=0.
**There is NO special area/volume correction, NO L'Hôpital series form, NO
sub-cell balance, NO reflective ∂ψ/∂r=0 stencil for the singular cell.**
The ONLY special r=0 treatment is the μ=−1 starting-direction sweep
(3.432)-(3.435) — and that is an ANGULAR-cascade SEED (gets the central
flux VALUE right), not a spatial-order remedy. So standard textbook
diamond/WDD spherical SN does NOT carry an O(h^2)-at-origin spatial
construction; the central cell inherits whatever order plain DD gives on
a cell whose inner face has zero area — which is degraded (first-order is
consistent with the well-known DD-on-singular-geometry behaviour).

## Bailey-Morel-Chang (2010) NSE 165:149-169 — LOCAL, ANGULAR only

`scratch/literature/Bailey-Morel-Chang(2010)...pdf`. VERIFIED by PDF read.
DOI 10.13182/NSE08-66 (NOT NSE08-64 — fix any stale cite). Their analysis
keeps r CONTINUOUS (Eq. (10) is angularly-discretized only); they
explicitly say (p.153, after Eq. 32): "we will not see an error in the
flux at the origin ... when the system is sufficiently thick and diffusive"
— this "error at the origin" is the angular FLUX DIP, not a spatial
truncation order. The paper CANNOT speak to spatial r=0 order. The M-M
weighted-diamond and its "first-order diffusion-limit preservation" is the
ANGULAR redistribution fix, orthogonal to spatial mesh order.

## Lathrop (2000) NSE 134(3):239-264 — the brief's "decisive" paper, NOT spatial

**Full ref**: K.D. Lathrop, "A Comparison of Angular Difference Schemes for
One-Dimensional Spherical Geometry S_N Equations," NSE 134(3):239-264 (2000).
**DOI 10.13182/NSE00-A2114** (NOT A2113). OSTI biblio/20015675. NOT in
scratch/literature/ — PAYWALLED (T&F); OSTI is metadata-only.
⚠ The brief framed this as a spatial-order paper — IT IS NOT. OSTI abstract
verbatim: *"To isolate errors caused by angular differencing, the
approximations are derived from the transport equation WITHOUT spatial
differencing, and the resulting coupled ODEs are solved with an ODE
solver."* It compares 5 ANGULAR schemes (diamond, M-M weighted diamond,
linear-continuous = original S_N, linear-discontinuous, new
quadratic-continuous). Its central-cell remark: *"the desirability is shown
of using an initializing computation of the μ=−1 angular flux to correctly
compute the central flux"* — this is the SEED-VALUE (starting-direction)
correctness point = Hébert (3.432)-(3.435), an angular issue. **Acquiring
Lathrop 2000 will NOT answer the spatial-order question** — it deliberately
zeroes spatial error. Do not chase it as the decisive spatial source.

## Wu, Zhongsheng & Fischer (1999) NSE 133 — the actual spatial remedy lead

"A Discrete Ordinates Nodal Method for One-Dimensional Neutron Transport
Calculation in Curvilinear Geometries," DOI 10.13182/NSE99-A2095, NSE 1999,
68 citations, PAYWALLED (not local). Abstract: Green's-function nodal SN
with Legendre spatial expansion; *"very high precision on coarse spatial
meshes relative to the standard fine-mesh SN method with the spatial
diamond-differencing scheme."* This is the class of method that beats DD's
spatial order in 1-D curvilinear — a nodal/higher-order or
characteristic/step-characteristic spatial scheme is the documented route
to better-than-DD central accuracy. Best single acquisition target if a
published O(h^2)-at-origin spatial closure is wanted.

## astro-ph starter-intensity papers — confirm r=0 is a SEED problem

Caro Gómez (UPM) 2002, arXiv astro-ph/0207291 "Inverse square law for the
starter intensity of the discrete ordinates method in spherical geometry"
(saved locally via webfetch). Companion astro-ph/0207652. These treat the
central flux as a STARTER-INTENSITY / starting-direction consistency
problem (μ=−1 transfer eqn → central flux obeys inverse-square law; isotropy
condition at center). Independent confirmation that getting the central
flux VALUE right is an angular-seed issue, not a spatial-order fix.

## Production-code reality

Hébert §3.9.4 IS the DRAGON formulation (uniform DD + μ=−1 seed, no singular
-cell spatial correction). PARTISN/DANTSYS lineage (Lathrop-Carlson DTF/
TWOTRAN) use diamond/weighted-diamond spatially with set-to-zero or
fixup negative-flux handling — no published special O(h^2) central-cell
spatial stencil. Conclusion: first-order-at-origin under plain DD is the
ACCEPTED de-facto behaviour in production curvilinear SN; the documented
way to recover high central accuracy is a higher-order SPATIAL scheme
(nodal / linear-discontinuous / characteristic), NOT a patched DD cell.

## Bottom line for ORPHEUS

- NO canonical reference reports a drop-in O(h^2) diamond/WDD central-cell
  SPATIAL closure. The first-order-at-origin is real, under-documented as a
  spatial-order statement, and consistent with DD on a zero-inner-area cell.
- The μ=−1 starting-direction sweep (Hébert 3.432-3.435) fixes the central
  flux VALUE/dip (angular), and ORPHEUS already targets it (Phase D) — it is
  NOT the spatial-order lever.
- To recover O(h^2) AT r=0, the literature points to a higher-order SPATIAL
  scheme on the singular cell: step-characteristic / linear-discontinuous /
  nodal (Wu 1999), or a sub-cell balance. None is a "corrected pole-face
  area" within plain DD.
