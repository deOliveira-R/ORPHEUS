---
name: alpha-dome-chart-vs-measure-cross-method
description: P4 — the SN α-dome is the discrete image of a CHART object (Pomraning 1989 Eq. 68, the second fundamental form of the sweep level surfaces); the quadrature is only the reconstruction rule. Per-method consumption table for SN/Pn/SPn/diffusion/MoC/CP/MC + DG-in-angle/Fn/kinetics/adjoint/FP/DSA, each with its structural reason.
metadata:
  type: project
---

# α-dome: chart, measure, or discretisation — and who could ever consume it

Attack 2026-08-27, deliverable `scratch/p4_alpha_cross_method.md`. Read-only.
**Why:** an architecture decision hung on which solver families could consume
`alpha_dome`. **How to apply:** fire on ANY "should this primitive be shared /
where does this geometric coefficient live" question in curvilinear transport,
and on any cross-method "who consumes X" brief.

## ⭐⭐⭐ THE ONE EQUATION THAT DECIDES EVERYTHING — Pomraning 1989 Eq. (68)

`scratch/literature_ocr/Pomraning(1989)The Transport Equation in General Geometry.md`
(NSE 101:330). His Eq. (63) says the angular-derivative terms ARE the gap between
`d/ds` and `∂` at fixed local angles; Eq. (68) gives the coefficient for ANY chart:

```
dμ/ds = −(1−μ²)(cos²φ/ρ_u + sin²φ/ρ_v)  −  [μ(1−μ²)^½/h_w]·[(cosφ/h_u)∂_û h_w + (sinφ/h_v)∂_v̂ h_w]
```

`ρ_u, ρ_v` = **principal radii of curvature** of the sweep coordinate's level
surface; `h_*` = **metric coefficients**. The first term is `−II(Ω_⊥, Ω_⊥)`, the
**second fundamental form / shape operator**. ⟹ **ZERO quadrature, ZERO
discretisation on the RHS.** The redistribution coefficient is pure extrinsic
differential geometry.

⭐ It DERIVES three previously-asserted rulings from one place:
* **slab** `ρ_u=ρ_v=∞` ⟹ `dμ/ds = 0` — the neutral element is a theorem.
* **sphere UMBILIC** `ρ_u=ρ_v=r` ⟹ `cos²φ+sin²φ=1` ⟹ **φ cancels** ⟹ one dome, no
  level structure.
* **cylinder FLAT ALONG ẑ** `ρ_v=∞` ⟹ the only channel that could move `μ_z` is
  identically zero ⟹ **#235's axis pair is wrong** is now a vanishing-principal-
  curvature theorem, not an observation.

⚠ Pomraning's normal-orientation convention makes `ρ` signed — quote the
STRUCTURE, not his signs.

## The verdict, three layers (do not collapse them)

| layer | object | owner |
|---|---|---|
| continuous | `f(x) = r·dx/ds` = `−r·II(Ω_⊥,Ω_⊥)` + lapse-gradient | **the CHART** |
| target | `f(x_{m±½})` | chart **evaluated on** the PARTITION |
| shipped recursion | `α = −Σ_{k≤m} w_k x'_k` | **the QUADRATURE** — reconstruction only |

⟹ **α is a CHART object reconstructed through the MEASURE for a DISCRETISATION
reason** (cumulative sum ⟹ discrete balance TELESCOPES exactly; pointwise `f(e)`
would balance only to quadrature accuracy).

⭐ **First split "angular measure" into (a) the level's INVARIANT measure — chart
data — and (b) the QUADRATURE — a discretisation choice.** The brief's trichotomy
is undecidable until you do; `f` depends on (a), `alpha_dome` on (b).

Two tests that point OPPOSITE ways (why it is a verdict, not a preference):
hold measure fixed / change chart ⟹ α changes (not a measure invariant); hold
chart / refine measure ⟹ α(e) has a **continuum limit as a function of the angle**
(so it IS a chart object).

## `f = r·dx/ds` per arm, and the ½ nobody had derived

`x` = the coordinate in which the LEVEL's own invariant measure is flat.

| | sphere | cylinder |
|---|---|---|
| `x`, measure | `μ`, `dΩ=dμ dφ` | `ω`, `dΩ=dμ_z dω` |
| `f` | `1−μ²` | `ξ = sinθ_p sinω` |
| `df/dx` vs increment `−w x'` | `−2μ` vs `−μ` ⟹ **factor ½** | `+η` vs `−η` ⟹ **factor 1** |
| `α = κ·f(edge)` EXACT? | **no**, asymptotic (Gauss in the COSINE, edges by cumulative WEIGHT) | **YES at every M** |

⭐ **The ½ is `df/dμ = −2μ`** — it decodes the tree's undocumented
"factor-of-2-absorbed normalization" AND the sphere-only `scale = 0.5` at
`angular_differencing.py:326`.

`[D]` **Cylinder closed form, exact:** nodes `ω_m=(2m+1)π/2M`, equal weights,
ω-midpoint edges `ω_{m+½}=(m+1)π/M`; Dirichlet sum
`Σcos((2k+1)θ)=sin(2(m+1)θ)/2sinθ` gives `α = −C_M·ξ(ω_{m+½})`,
`C_M = W/(2M sin(π/2M)) → W/π`. ⟹ **the dome closure IS `sin π = 0`** and the
"dome" is literally the graph of `sin`. A tolerance constant + guard + `ValueError`
exist for an identity.

⛔ **RETRACTED IN-ATTACK:** "both arms have the same flux function `sin(angle)`;
`1−μ²` vs `ξ` is a chart artefact." FALSE — the coordinates differ because the
RESIDUAL SYMMETRY GROUPS differ. The unification is the FORMULA `f = r·dx/ds`,
never a common closed form.

## Per-method table — consumed? / the structural reason

| family | ? | reason |
|---|---|---|
| **SN** | YES | the only value-path consumer; finite-volume in angle |
| **Pn** | no | `μ·` and `(1−μ²)∂_μ` are components of ONE rank-1 tensor ⟹ same `ℓ→ℓ±1` rule. `∂_μ[(1−μ²)P_ℓ] = [ℓ(ℓ−1)P_{ℓ−1} − (ℓ+1)(ℓ+2)P_{ℓ+1}]/(2ℓ+1)` folds into the recurrence Pn already has. `[L]` Garcia 2019 JCP **Eq. (3)** verbatim |
| **SPn** | no | DEFINED by `d²/dx² → ∇²` on the **PLANAR** Pn set — the plane has no term to carry. Geometry moves to the SPATIAL base. `[L]` Makine 2018 Eqs. 13–15 (all `∇²`) |
| **diffusion** | no | `ℓ(ℓ−1)` **vanishes at ℓ=0,1**. ℓ=0 row ≡ 0; ℓ=1 reaches only φ₂, killed by P1. What survives is the `2φ₁/r` that IS the spherical divergence ⟹ Fick has **no 1/r term** |
| **MoC** | no | Ω constant in the **GLOBAL** frame ⟹ `Ω·∇ = d/ds`, chart-free. Never generated. Relocates to the TRACK SEGMENTATION table. `[M]` `orpheus/moc/` = **0** `CoordSystem` hits |
| **CP** | no | angle integrated out analytically ⟹ no angular unknown. **α is the cost of COLLOCATING; CP quadratures the resolvent instead** |
| **MC** | no | Cartesian Ω, constant between collisions. No local frame ⟹ no connection |
| **DG/FE in ANGLE** | **YES — the genuine 2nd consumer** | needs `f` **pointwise at ITS edges + in the interior**; needs no recursion |
| **Fn / Case** | no | `Φ = rΨ` makes it **formally PLANE geometry** with the term demoted to a SOURCE + Volterra. `[L]` Garcia 2019 Eqs. 9–11. `[M]` ORPHEUS's own sphere refs already do this (slab algebra + parity flip `s=∓1`) |
| **nodal** | inherits | spatial axis; `R = R_spatial ⊗ A_angular` is a MEASURED factorisation |
| **kinetics/α-eig** | no | `+α/v` is a multiplication operator. ⚠ **one-letter collision** — discriminate by INDEX SET + CO-OPERATOR |
| **adjoint** | yes | same operator transposed; not a new consumer |
| **Fokker–Planck/CSD** | consumes `f` | `∂_μ[(1−μ²)∂_μ]`, identical `f`, one order up |
| **consistent DSA** | consumes α as DATA | Schur complement of the ASSEMBLED high-order op |

⭐ **The 3-clause discriminator that generates every row:** a family needs α iff it
has (1) an angular UNKNOWN, (2) indexed by a LOCAL rotating frame, (3) whose
derivative is discretised by COLLOCATION. Every "no" fails a nameable clause.

⭐ **ONE theorem, three representations** (`f(±1)=0` ≡ Pn's `ℓ(ℓ−1)|_{ℓ≤1}=0` ≡
`α_{1/2}=α_{M+½}=0`). Use it to convert "diffusion approximates it away" into
"the operator ANNIHILATES the retained subspace" — an exact algebraic claim.

## The generalisation verdict — and the L60 trap, on a FUNCTION

**The shared primitive is `f`, NOT `alpha_dome`.** `alpha_dome` = `f` ∘ three
invisible discretisation decisions: collocate-in-angle / these EDGES / reconstruct
by cumulative quadrature of `f'`.

L-004's test: `f` has ≥2 non-iso consumers with **3 genuinely different applied
morphisms** (cumulative-quadrature | pointwise | `⟨P_k, ∂_μ(fP_ℓ)⟩`) ⟹ shared.
`alpha_dome` has ONE consumer and the only morphism is `id` ⟹ SN-specific.

⭐⭐ **L60 on a function signature.** `alpha_dome(ndarray, ndarray)` has free-name
set `{np}` — maximal apparent generality. **But the arrays ARE a quadrature**, and
every non-consumer lacks that object by construction. The signature records the
author's RECONSTRUCTION CHOICE and reads as a claim about the object's nature.
**Discriminating NEGATIVE test:** write what a DG-in-angle caller needs — `f(e)` at
ITS edges and `f(x)` in the interior. `alpha_dome` supplies **neither** (no edge
argument; output defined only at the SN partition). **A second consumer cannot
call it.** REDs today; cannot be passed by an implementation I call correct.

## Byproduct — the `k` the ladder names and the diagnostic drops

`[M]` P4 (`angular_differencing.py:91`) is `α_{m+½} = α_{m−½} − k·w_m·η_m` with an
explicit `k`; production hard-codes `k=1` (`reduced_operator.py:842`);
`alpha_defect_beta` (`:419`) grades that against Lathrop's `k=2` target `(1−e²)`.
`[D]` with `Σw=2` (pinned by `diffusion_limit_c`'s own `[M]` `Σwμ²=0.666666666667`)
the `k=1` dome → `(1−e²)/2`, so the "defect" is dominated by `−(1−e²)/2`, O(1).
`S₂` GL: `α=(0,0.5774,0)`; vs `f/2` → `+0.077`; vs `f` → `−0.423`. Analysis-only
(L0), no production value moves. NOT acted on (read-only brief).

## Refuted here, with reasons (do not re-attack)

Riemannian geodesic flow (ambient space is FLAT — the object is EXTRINSIC curvature
of the level surfaces, so "Christoffel symbols of space" is a category error) ·
optimal transport (signed linear operator, no cost functional, no marginal) ·
angular homogenisation (no small angular parameter; the thick limit is SPATIAL) ·
Krylov on the cascade (triangular, one pass, no iteration) · Magnus/Lie-group (a
scalar two-term relation, no operator ordering) · sphere→slab reduction FOR SN
(exact, but trades a triangular sweep for a Volterra inner iteration — wrong
direction for a sweep solver; it is the explanation of the Fn row, not a proposal).
Plus the standing set: DEC/homology, MPO (bond dim 1), Wiener–Hopf, category
theory, Miller/Gautschi (`G(M)=1`), MaxEnt.

## Open

Morel (1985) *An improved Fokker–Planck angular differencing scheme* — **NOT in
`scratch/literature/`** (the folder has Morel 1982 DSA and Morel 1989 hybrid
collocation–Galerkin). Acquisition must be REQUESTED, not assumed. It is the
published precedent for the ladder on the sibling operator, and reproducing
ORPHEUS's α from it would satisfy §5.2's ≥2 criterion by evidence.
