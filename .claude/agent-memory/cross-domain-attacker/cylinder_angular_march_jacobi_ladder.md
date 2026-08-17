---
name: cylinder-angular-march-jacobi-ladder
description: Q68/#235 — the curvilinear angular redistribution IS the first-order Jacobi–Sturm–Liouville factor of the level's own generating measure; the level's ordinates ARE that measure's Gauss nodes (cylinder = Chebyshev–Gauss, sphere = Gauss–Legendre) ⟹ tridiagonal ladder, free truncation, no seed. Plus: the operator admits NO boundary condition.
metadata:
  type: project
---

# The cylinder's angular march — the native frame, and why the seed is an artifact

Attack on "is the 1-D η-march the right structure for the cylinder's angular
variable?" Deliverable `scratch/q68_angular_frame_attack.md`. Branch-verified
reads (L-005) on `refactor/operator-strategy-layers`; no shell, so every number
is either an in-tree `[M]` citation or symbolically derived `[D]`.

**Why:** #235 is filed as "cylindrical anisotropic accuracy needs a 2-D (η,φ)
angular closure". The frame attack says the axis pair is wrong AND names the
native 1-D object. **How to apply:** fire on any curvilinear angular-differencing,
seed, α-dome or τ question in either geometry; the first three rulings generalise
past SN.

## ⭐⭐⭐ THE NATIVE FRAME — one sentence, both geometries

> The curvilinear angular-redistribution operator is the **first-order factor of
> the Jacobi–Sturm–Liouville operator of the level's own generating measure**,
> and the level's ordinates are **that measure's Gauss nodes**.

`[D]` The conservative redistribution is `∂_x[(1−x²)w(x)·]` with `w` the level's
Gaussian weight — Legendre (`w=1 ⟹ 1−x²`) on the sphere, Chebyshev-1
(`w=(1−x²)^{−1/2} ⟹ (1−x²)^{1/2}=ξ`) on a folded cylinder level. That is
`jacobi(0,0)` and `jacobi(−½,−½)`, both already in
`orpheus/numerics/generating_measure.py`.

`[D]` **The shipped cylinder nodes ARE Chebyshev–Gauss.** `folded_product` =
GL(n_μ) × trapezoid(n_φ, **STAGGERED**) folded by σ_y ⟹ `ω_k = (k+½)π/M`,
`M = n_φ/2` ⟹ `η_k/sinθ = cos((2k+1)π/(2M))` = the roots of `T_M`, equal
weights = Chebyshev–Gauss weights, degree `2M−1`. **1 of 1** admitted cylinder
rule families (`assert_carrying_quadrature` refuses the rest).

`[D]` **The operator is TRIDIAGONAL with zero diagonal** in the family's basis:
* cylinder: `A cos(kω) = ½(k+1)cos((k+1)ω) − ½(k−1)cos((k−1)ω)`
* sphere: `∂_μ[(1−μ²)P_ℓ] = (ℓ(ℓ−1)P_{ℓ−1} − (ℓ+1)(ℓ+2)P_{ℓ+1})/(2ℓ+1)`

Three theorems fall out, all `[D]`, all one line:
1. **Exact particle balance at EVERY truncation** — the `j=0` row is identically
   zero (the down-step coefficient vanishes at `k=1`/`ℓ=1`). Lathrop's stated
   design requirement, obtained structurally instead of by α-recursion + a
   runtime closure guard.
2. **The truncation is FREE** — the top mode's up-step produces `T_M`/`P_N`,
   which **vanishes at its own Gauss nodes**. Zero collocation error, no ad-hoc
   top-mode closure.
3. **The transpose is `A_ord.T`** — the hand-rolled reverse recurrence dies.

Deleted, not moved: `angular_cell_edges_per_level`, `morel_montry_tau_per_level`,
`alpha_dome` + `_assert_alpha_dome_closes` + `_ALPHA_CLOSURE_ATOL`, `MarchStart`,
`march_start_structure_per_level`, `assert_carrying_quadrature`,
`edge_extrapolated_seed`, `_psi_half_grid_single_level`, the `angular_adjoint`
recurrence.

⛔ **Three honest costs.** (a) The sweep loses angular triangularity — per-cell
`M×M` solve with a precomputable factorisation; fine at `nx≤320, M≤32`, **does
NOT transplant to 2-D/3-D**. (b) Gibbs on ω-nonsmooth ψ (annulus tangent rays,
directional sources) — land it as another `PoleAngularClosureBase` registry
member, not a replacement. (c) ⛔ **the published failure mode applies**:
Larsen–Morel 2010 §1.5.1 — LD-in-angle LOST to weighted diamond precisely
because it *"decouples from the starting-direction flux"*, which is computed
more accurately than the rest. A seedless spectral scheme discards both accurate
legs. Fix = the constrained-Galerkin synthesis below.

## The second structural finding — the operator admits NO boundary condition

`[D]` `∫₀^π ∂_ω(ξψ)φ dω = [ξψφ]₀^π − ∫ ξψφ' dω`, and **the boundary term
vanishes identically** because `ξ(0)=ξ(π)=0`. Fichera function zero at both
ends: both endpoints are characteristic; neither is inflow or outflow. This is a
degenerate first-order operator (Fichera / Oleĭnik–Radkevič), **not an IVP**.

⟹ **the seed `ψ_{1/2}` is a discretisation artifact, not physical data.** What
the endpoints supply is a **compatibility condition**, and there are TWO
(the plain radial ODE at `η=∓sinθ`) — exactly the two legs
`RadialCharacteristicOperator.solve` already marches. `[D]` the discrete count:
`2M+1` unknowns vs `2M` equations ⟹ **the scheme is ONE CONDITION SHORT and TWO
DATA LONG**, while the continuous problem is well-posed with none.

⭐ **The synthesis that resolves the seedless-scheme risk**: expand in `{T_k}`,
impose the two marched leg values as endpoint constraints
(`Σ(−1)^k ĉ_k = ψ(π)`, `Σ ĉ_k = ψ(0)`), drop the two top balance equations.
That is **Lanczos' tau method**, it USES both accurate legs, and it is exactly
Z₂-symmetric because the two constraints are each other's reflection images.
⚠ **NAME HAZARD:** `tau` is already overloaded 3× in-tree; spell it
`endpoint_constrained`, never `tau`.

## Rulings that are theorems, not opinions

- **α IS the flux function.** `[D]` `α(ω) = (W/π)·ξ(ω)` — reproduces the memo's
  measured `α ∝ ξ(e_arc)` to 1e-16. The "dome" is the graph of `sin ω`;
  `α_{1/2}=α_{M+1/2}=0` is `sin 0 = sin π = 0`. A tolerance constant, a guard and
  a raise exist for an identity.
- **⛔ #235's axis pair is WRONG.** `Ω·∇ = η∂_r − (ξ/r)∂_ω` has **no `∂/∂μ_z`
  term** — `μ_z` is exactly conserved (ẑ-translation), levels are exactly
  decoupled, the redistribution is exactly 1-D. There is no angle-angle stencil
  to build; the literature survey's "found NOTHING published" is consistent with
  "there is nothing to publish". **Two legitimate two-variable readings survive**:
  the CHART question (τ barycentric in η, cells+nodes uniform in ω) and the
  space-angle `(r,ω)` coupling along `p = r sinω`.
- **One parameter, two conditions.** τ must be both the β=0 barycentre (BMC-43)
  and `½+O(w)` (R&L 15/16). Lathrop 2000 Eq. (30) is the published statement
  that it cannot be both. **The cell is one DOF short** — that is a countable
  claim, and it is what Lathrop's schemes 4/5 (LD-in-angle, quadratic-continuous)
  and the published `quadratic-first-cell + LD-rest` hybrid actually fix.
- **The `Z₂` the method breaks.** `R: ω↦π−ω` gives `R A R = −A` exactly; the
  closure family is anti-equivariant, the one-sided march is not. `R` acts on
  `cos mω` with `(−1)^m`, so the odd sector's leading member is `cos ω ∝ η` =
  **the radial current** — the moment an anisotropic MMS sees.
- **The existing gate is blind and the operator identity is the fix.** The tree
  gates `Π(1−τ)/τ = 1`, a SCALAR — the `(0,M)` entry of the operator statement
  `R A_h R = −A_h`. `[M-cite]` all four τ variants (shipped/diamond/reversed/
  shuffled) pass the scalar gate. ⚠ the operator probe only separates them at
  **M ≥ 4** (at M=2 reversal preserves `τ₁+τ₂=1`).
- **Smell #16 Shape 1 at maximum confidence.** `outgoing_face_from_average` and
  its VJP (`transport/spatial/scheme.py:1197, 1215`) are **character-identical,
  same reduction tree**, to the angular recurrence and its adjoint
  (`pole_angular_closure.py:1184, 1676`) — and
  `psi_half_angle_seed.py:105` **already imports the shared VJP for the radial
  march**. One module family imports the primitive for one march and
  re-implements it 1500 lines away for the other. Collapsing it makes
  `LinearDiscontinuous`-in-angle (Lathrop scheme 4) *reachable* — it is already a
  registered member of the spatial family, never instantiated on the angular axis.

## The instrument the campaign is missing

`[M-cite]` `q64_tau_partition_memo.md` §0a: 3 of 5 instruments dead, the
trajectory-resolvent cross-check **reference-limited at ≈3e-2**. The
reference-free replacement is one closed form `[D]`: a **purely-absorbing
cylinder**, `ψ = (Q/4πΣ_t)(1−exp(−Σ_t L/sinθ_p))` with
`L = r cos ω_m + √(R² − r² sin²ω_m)` (checks: `ω=π ⟹ L=R−r`; `ω=0 ⟹ L=R+r`).
No quadrature, no iteration, no scattering, no reference solver. Refine `nx`
until the residual plateaus and read pure angular-closure error per ordinate.
This is Lathrop 2000's own construction, transplanted to the cylinder.

## Test order (what to run, in what order, and why)

1. **Free, today** — `R_p = ψ_{M+1/2}^{marched} − cells(p,+1)`. Both arrays
   exist. Predict order 1 and correlation with the MMS error; the τ≡½ mutation
   must move both together. Falsifier: order 2, or anti-correlation.
2. **Decides #235** — order-FIT (not ratio) of the aniso MMS in `n_φ` at
   `nx=320` for M-M τ vs τ≡½. Hypothesis A (Lathrop-30, δ≠0 ⟹ first order):
   1 vs 2. Hypothesis B (genuine 2-D gap): both flatline. Mutually exclusive.
   ⚠ the campaign's measured **2.5–3.0× ratio is consistent with A but is not
   evidence** — a ratio is not an order.
3. **Bit-identity** — single-source the closure onto the spatial primitives;
   `array_equal`, 0 ULP. A RED means `angular_adjoint` is not the exact
   transpose and the #208 dagger claim over the angular block is weaker than
   advertised.
4. **P5 exactness** — extend the `angular_differencing.py` P0–P4 ladder with
   "exact on `cos 2ω`". Both schemes pass `k=0,1`; only the ladder passes `k=2`.

## REFUTED with reason (do not re-attack)

- **Miller/Gautschi minimal-solution ("march the stable direction")** — ⛔ `[D]`
  `G(M)=1` exactly ⟹ the recurrence is **reversible**; both directions carry the
  identical amplification profile `G(m)`. No dominant/minimal pair, no purchase.
  **High-prior trap**: the measured `A(M)=9.44` actively baits it.
- **Two-point BVP** — superseded, not refuted: there are no boundary conditions,
  so "BVP" is the wrong noun. Constrained Galerkin, two *compatibility* conditions.
- **Step-characteristic / exponential in angle** — a **measured** dead end
  (Lathrop: `O(Δμ)`, decouples from the seed, hybrid shows first-order
  absorption). And the α-weighting does NOT secretly encode it: the ω operator is
  purely conservative, with **no zeroth-order sink**, so there is no exponential
  attenuation in ω to encode.
- **MPO / tensor networks** — bond dimension **1** (`[M]` the recurrence never
  mixes cells).
- **DEC / de Rham** — 1-D complex ⟹ `H¹=0`; the only content is telescoping,
  already had.
- **Symplectic** — the conserved quantity IS `p = r sinω`, and the chart that
  conserves it by construction IS the characteristic chart = MoC. No independent
  lever.
- **DSA to carry the diffusion limit** — ⛔ DSA changes the ITERATION, not the
  FIXED POINT; `β≠0` is a property of the converged solution.
- **MaxEnt** — the shape is fully resolved by M ordinates; the closure is a
  reconstruction between resolved nodes, not a moment closure. No entropy deficit.

## Pollination ACQUIRE — the highest-value open lead

**Morel's Fokker–Planck angular differencing** (Morel 1985, NSE, "An improved
Fokker-Planck angular differencing scheme"). The FP operator
`∂_μ[(1−μ²)∂_μ·]` has the **same degenerate endpoints and the same
Sturm–Liouville factor**, one differential order up, and Morel's scheme is built
to be exact on the Legendre eigenbasis — i.e. the published precedent for the
ladder construction, on the sibling operator.
`scratch/q68_curvilinear_angular_differencing_survey.md` flags it and defers to a
**§Q1.5 that was never written** (the file ends at §Q1.4). One web probe found
only *spatial* Chebyshev-collocation RTE work; the angular-side construction is
**not placed** — do not claim novelty.
