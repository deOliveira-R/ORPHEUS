---
name: Sanchez-McCormick §III.F.1 rank-N per-face interface currents
description: Canonical formulation of multifunction per-face interface currents for CP — moment basis, W matrix, P_Si escape, A_αβ albedo, closure, conservation, reciprocity — for comparison against ORPHEUS rank-N hollow-sphere primitives.
type: reference
---

# Sanchez & McCormick 1982 §III.F.1 — canonical rank-N per-face interface-current formulation

Source: R. Sanchez & N. J. McCormick, *Nucl. Sci. Eng.* **80**, 481-535 (1982). PDF at `/workspaces/ORPHEUS/1982NSE80-481.pdf`. Section §III.F is "Nodal Methods" (p. 520); §III.F.1 "Interface Current Method" occupies **pp. 520-523, equations (162)-(169)**. The paper's only explicit discussion of rank-N per-face interface currents is this block. The transmission/albedo operator is introduced via Eq. (10) earlier in §I; the expansion philosophy (multifunction CP, shifted Legendre moments) is inherited from §III.C.4 / §III.E. Table II (referenced for the angular factor ω) lists the three 1D geometries with α = 0 (slab), 1 (cylinder), 2 (sphere).

## 1. Section location — confirmed

- **§III.F Nodal Methods** starts mid-p. 520 (right column).
- **§III.F.1 Interface Current Method** — p. 520 right col through p. 523 top. Principal equations (162), (163), (164), (165), (166), (167), (168), (169).
- §III.F.2 "Transverse Nodal Method" follows on p. 523 and is *not* relevant here.

## 2. Moment basis — definition (Eq. 165)

The angular flux on each subsurface α of cell is expanded in **surface modes** f^ρ_{±,α}(r_b, Ω):

```
ψ_±(r_b, Ω) = Σ_{α,ρ}  J^ρ_{±,α} · f^ρ_{±,α}(r_b, Ω)                        (165)
```

ρ = 1…P_α indexes the mode within subsurface α. f^ρ_{±,α} is defined on subsurface α for ±Ω·n > 0 and is zero otherwise (half-range basis — one basis per emission half-space per subsurface).

The normalization is the weighted-orthogonality condition (unnumbered inline, between 165 and 166):

```
∫ f^ρ_{+,α}(r_b, Ω) · f^ν_{+,α}(r_b, Ω) · (Ω·n) dΩ dA  =  (π A_α)^{-1} · δ_αβ δ_ρν
```

where A_α is the area of subsurface α. The projection inner product carries a **µ-weight (Ω·n)**, so the moments J^ρ_{±,α} are **partial-current moments** (µ-weighted), **not plain angular-flux moments**.

**First mode ρ = 0 is taken constant:** f^0_{±,α} = (π A_α)^{-1}. Therefore (footnote-style sentence on p. 521):

> *"the first expansion function is taken to be a constant, which gives f^0_{±,α} = (π A_α)^{-1}; hence, J^0_{±,α} is the total neutron current leaving (+) and entering (-) the cell through subsurface α. For this reason, the J^0_{±,α} may be called the IC components."*

This is the load-bearing normalization: **J^0 = ∫_{Ω·n>0} Ω·n ψ dΩ dA** — the total outgoing current *per subsurface*, **not** the current density, **not** divided by A. Higher modes are the corresponding higher µ-weighted angular-flux moments.

Gelbard (2n+1) placement: Sanchez-McCormick do **not** put (2n+1) in the basis functions themselves. The mode basis is orthonormal under the µ-weighted inner product with Kronecker delta (no 2n+1 in the RHS of the orthogonality relation). Any (2n+1) factor therefore appears *only* when Sanchez-McCormick use the shifted-Legendre form P̃_n as the spatial f_α; the basis weight itself is (π A_α)^{-1} × δ.

## 3. Transmission / surface-to-surface matrix — Eq. 166 cluster

The three defining surface coefficients in the boxed equation set following (166):

```
P^{kρ}_{iS_α} = -4π V_i^{-1} ∫ f^k_i(r) dr ∫ g(r'_b → r) · f^ρ_{-,α}(r'_b, Ω_s) · (Ω_s·n) dA'   [volume ← incoming surface, G_bc analog]

P^{ρk}_{S_α i} = +π A_α · V_i^{-1} ∫ f^k_i(r) dr ∫ g(r → r'_b) · f^ρ_{+,α}(r'_b, Ω_s) · (Ω_s·n) dA'  [P_esc analog — volume → outgoing surface]

P^{ρν}_{S_α S_β} = -4π² A_α ∫ f^ρ_{+,α}(r_b, Ω_s)(Ω_s·n) dA
                          × ∫ g(r'_b → r_b) · f^ν_{-,β}(r'_b, Ω_s)(Ω_s·n') dA'   [W transmission — surface → surface]
```

The collision kernel is `g(r'→r) = exp(-τ) / (4π s²)` (Eq. 107), and angular streaming appears through the arguments of g evaluated at line-of-sight.

Normalization summary per primitive:

| Primitive | Leading scalar | Source volume | Arrival normalization |
|---|---|---|---|
| `P^{kρ}_{iS_α}` (G_bc analog — incoming surface → volume mode k) | −4π V_i^{-1} | volume V_i **divides** | dA' weighted by (Ω·n) |
| `P^{ρk}_{S_α i}` (P_esc analog — volume mode k → outgoing surface ρ) | +π A_α · V_i^{-1} | V_i divides | dA weighted by (Ω·n) |
| `P^{ρν}_{S_α S_β}` (W transmission — incoming β ν → outgoing α ρ) | −4π² A_α | (none — pure surface-surface) | dA on emission A_α, dA' on incidence A_β, both with (Ω·n) |

**Critical observation 1:** the surface-to-surface W carries a **leading −4π² A_α** *on the emitting subsurface*. It is **normalized per unit of cell-total outgoing current**, not per unit area. The sign is absorbed by the convention (Ω_s · n) with n outward on the emission face.

**Critical observation 2:** there is **no** `(ρ_max/R)²` Jacobian in these Sanchez-McCormick integrands. The surface integrals are performed directly in (r_b, dA, Ω) — streaming geometry is encoded in the collision kernel g(r'_b → r_b) through `s = |r_b − r'_b|`, not in a separate surface Jacobian. The `(ρ/R)²` factor is ORPHEUS's own addition to map observer-centered (r_i, Ω) integration to surface area, and it belongs *only* in the single-surface `compute_P_esc_mode` form (which uses observer-centered coordinates). Once the per-face decomposition is done with *surface-centered* coordinates, the Jacobian conflates with the surface element dA and must not be re-applied.

## 4. Per-face escape / G_bc analog

Already embedded in the Eq. (166) cluster above:

- **`P_esc` analog** is `P^{ρk}_{S_α i}` — volumetric source → outgoing surface ρ at subsurface α.
- **`G_bc` analog** is `P^{kρ}_{iS_α}` — incoming surface ν at subsurface β → volume mode k at zone i.

Both carry (Ω_s · n) surface weight. Both use the free-space collision kernel g; no Jacobian. The leading factor difference ( −4π V^{-1} for G_bc vs +π A · V^{-1} for P_esc) is **4× and comes with the (Ω·n) weight**. It is not an artefact of normalization: it reflects that the incoming surface flux projects through 4π (full half-range integrated angular distribution), whereas the outgoing surface only sees one half-space (π).

**This exactly matches the F.4 (rank-1) scalar code** where `compute_G_bc_outer = 4 · compute_P_esc_outer` empirically — that's the ratio −4π V^{-1} / (π A V^{-1}) × A/A = 4.

## 5. Closure equation — Eq. (166) and reciprocity Eq. (167)

The full closed system for one cell is:

```
φ^k_i       = Σ_{j,l} P^{kl}_{ij} V_j F^l_j  +  Σ_{α,ρ} P^{kρ}_{iS_α} J^ρ_{-,α}                  (166a)

J^ρ_{+,α}   = Σ_{i,k}  P^{ρk}_{S_α i} V_i F^k_i  +  Σ_{β,ν} P^{ρν}_{S_α S_β} J^ν_{-,β}          (166b)

J^ρ_{-,α}   = J^ρ_{0,α}  +  Σ_{β,ν̄} A^{ρν}_{αβ̄} J^ν_{+,β̄}                                      (166c)
```

where Eq. (166c) is the boundary/interface coupling. For **white BC** (isotropic reflection across the cell boundary), the albedo operator simplifies and (164) gives J_0 = 0, reducing (166c) to J^-_α = A_αα J^+_α for each subsurface (diagonal in α because the reflection is the subsurface's own outgoing current).

**Closure form** (for an isolated cell with white BC):

```
(I − W · A) J^+  =  P_S · V · F                                   [matrix form]
```

where:
- W is the block matrix with entries P^{ρν}_{S_α S_β} (rows = outgoing (α,ρ), cols = incoming (β,ν)).
- A is the albedo block (diagonal in α for white BC on non-coupled subsurfaces).
- P_S is the block of P^{ρk}_{S_α i} (volume source → outgoing surface).

After solving for J^+, the volumetric flux is recovered by back-substituting J^- = A J^+ into (166a):

```
φ^k_i  =  Σ P^{kl}_{ij} V_j F^l_j  +  Σ P^{kρ}_{iS_α} (A J^+)^ρ_α                [166a with closure]
```

**Where (2n+1) appears:** Sanchez-McCormick (167) gives the reciprocity relations

```
P^{kl}_{ij}        = P^{lk}_{ji}                                    (167a)
P^{kρ}_{iS_α}       = 4 · (A_α)^{-1} · V_i · P^{ρk}_{S_α i}          (167b)   -- scalar part
A_β · P^{ρν}_{S_α S_β}  =  A_α · P^{νρ}_{S_β S_α}                   (167c)   -- NOTE TRANSPOSED MODE INDICES
```

Equation (167b) is the "4×" identity between G_bc and P_esc of §4 above. **(167c) is the generalized reciprocity — the mode indices ρ, ν transpose under the surface-swap.** Non-trivial for ρ ≠ ν and a test the W primitives must pass.

Gelbard (2n+1) factors only appear if one chooses shifted Legendre P̃_n as the spatial basis **and** projects under the unweighted inner product ∫ P̃_n P̃_m dµ = (2n+1)^{-1}. Sanchez-McCormick's orthogonality condition above is **µ-weighted with δ** (not δ/(2n+1)) — so the (2n+1) does *not* appear in their W or P_S. If a code adopts P̃_n directly with µ-weight already baked in (which is the Marshak DP_{N-1} convention), (2n+1) disappears. If a code computes P̃_n moments without the µ-weight (plain half-range Legendre), (2n+1) must then be reinserted in the reflection operator to recover the same physics.

**So the (2n+1) lives *exclusively in the reflection operator* when the basis is unweighted half-range Legendre, and *vanishes* when the basis is already µ-weighted half-range Legendre.** ORPHEUS's `reflection_marshak` carrying diag(1,3,5,…) is consistent with the unweighted-basis convention — but only *if* the P_esc / G_bc / W primitives are using the unweighted basis as well.

## 6. Differences from ORPHEUS's current per-face primitives

Ground truth ORPHEUS code at
`orpheus/derivations/peierls_geometry.py:1631-1790` (P_esc outer/inner modes) and 2921-3050 (W rank-N).

### 6a. The spurious (ρ_max/R)² Jacobian in `compute_P_esc_outer_mode`

ORPHEUS Eq. (docstring) as implemented:

```
P^{(n)}_{esc,out}(r_i) = C_d ∫ W_Ω · (ρ_out/R)² · P̃_n(µ_exit) · K_esc(τ) dΩ
```

But the actual code (line 1687-1712) **removed** the (ρ_max/R)² factor in the Phase F.5 convention update — comment says "Mark-Lambert angular measure (sin θ dθ) with P̃_n(µ_exit) Legendre factor — NO (ρ/R)² Jacobian". Good — this matches Sanchez-McCormick who have no such Jacobian in P^{ρk}_{S_α i}.

**BUT** the observer-centered integration is then:

```
P_esc^{(n)}(r_i)  ≈  pref · Σ_k  w_k · (sin θ_k) · P̃_n(µ_exit) · K_esc(τ_k)
```

which is the single-surface/Mark convention. Sanchez-McCormick's P^{ρk}_{S_α i} is instead **surface-centered**:

```
P^{ρk}_{S_α i}  =  +π A_α · V_i^{-1} · ∫_{V_i} dr ∫_{A_α} g(r → r'_b) · f^ρ_{+,α}(r'_b, Ω_s) · (Ω_s·n) dA'
```

i.e. integrate over subsurface points r'_b with a (Ω_s·n) weight — **that's the missing µ-weight on the surface basis**. ORPHEUS's observer-centered integral uses `sin θ dθ dφ` (Lambertian angular measure). Converting to surface-centered via Jacobian would give an extra `(Ω_s·n) = cos θ` factor that the code **has dropped along with the (ρ/R)² Jacobian**. That is the bug.

**The right Sanchez-McCormick form, observer-centered, carries `cos θ` (not `(ρ/R)²`) — they are not the same thing.** `(ρ/R)²` is the area Jacobian for observer→surface sphere-area projection. `cos θ = Ω·n` is the µ-weight in the surface basis inner product. The code dropped both, but only one is actually spurious — the other is essential.

### 6b. Model A vs Model B split

ORPHEUS's per-face N=1 code puts Model A (ray-blocking — exclude rays that hit inner) on P_esc and Model B (include cavity-crossers) on G_bc — and that reproduces F.4 bit-exactly. Sanchez-McCormick's §III.F.1 does not draw this distinction because they define the problem with *one* subsurface α per face and treat blocking entirely through the transmission matrix W = P^{ρν}_{S_α S_β}. In a hollow annular cell:

- Outer subsurface α = outer shell, f^ρ_{+,outer} = nonzero only on outer, Ω·n_outer > 0 points outward.
- Inner subsurface α = inner shell, f^ρ_{+,inner} = nonzero only on inner, Ω·n_inner > 0 points *into the cavity* (away from the cell medium).

A ray leaving the outer surface inward but hitting the inner surface first is captured by the `g(r'_b → r_b)` kernel in W between outer-outgoing and inner-incoming — it's neither P_esc nor G_bc, it's a W off-diagonal block. **So Sanchez-McCormick's formalism doesn't need a "Model A vs Model B" split at the primitive level** — the split only emerges when you use a *single* surface α = outer and try to collapse the inner into W in an ad-hoc way. ORPHEUS's N=1 scalar code does the collapse, which is why it must distinguish models.

This is a **structural mismatch**: ORPHEUS's current rank-N per-face code is trying to generalize the F.4 scalar (single-surface with manual model split) rather than implementing the native Sanchez-McCormick two-subsurface formulation. Every "Model A/B" split in rank-N is a Band-Aid over this structural choice.

### 6c. σ_t-dependent residual

From `peierls_rank_n_closure_not_basis.md`: optimal mode-1 scale `c_opt` grows with σ_t · R. That is the **fingerprint of a missing (Ω·n) = cos θ weight in the partial-current moment**. When σ_t · R is small (thin cell), nearly-grazing rays (cos θ ≈ 0) contribute heavily to the P̃_n(µ_exit) integral but should be down-weighted by cos θ for the partial-current moment. Without the cos θ, the mode-1 amplitude is systematically inflated; the inflation grows as σ_t · R makes the exponential attenuation milder (grazing rays get longer τ and cleaner exp(−τ) attenuation, making the missing cos θ bite harder when σ_t·R is small).

The scan found `c_opt ≈ 0.05 @ σ_t·R = 2.5`, `0.16 @ σ_t·R = 5`, `0.36 @ σ_t·R = 10`. The ratio `c_opt · (2n+1) = 1` is not achieved anywhere — confirming no constant fix. But `c_opt` varying with σ_t·R is consistent with a missing ∫ cos θ · P̃_n(µ) dµ weight (which integrates down exactly one of the (2n+1) factors — i.e. moving (2n+1) from the reflection operator into the basis µ-weight and replacing unweighted-basis with µ-weighted-basis).

## Bottom-line correction recipe

To align ORPHEUS's rank-N per-face primitives with Sanchez-McCormick §III.F.1:

1. **Drop the `(ρ_max/R)²` Jacobian everywhere** (already done in Phase F.5 for P_esc_outer_mode; confirm it's dropped in W and G_bc analogs too).
2. **Insert (Ω·n) = cos θ** into the **surface-centered** integration of every per-face primitive. In observer-centered form that becomes an extra `cos θ` (not `sin θ`) weight inside the dΩ integral for P̃_n modes with n ≥ 1. Mode 0 keeps the Lambertian form and remains bit-exact with F.4.
3. **Split by subsurface, not by Model A/B.** Define an outer-subsurface and an inner-subsurface basis; implement the six primitives (P_S_outer, P_S_inner, G_bc_outer, G_bc_inner, W_outer-outer, W_inner-outer, W_outer-inner, W_inner-inner — eight total). Blocking is a side effect of the surface-to-surface integral, not a model switch.
4. **Reflection operator stays diag(1, 2n+1)** under the µ-weighted basis interpretation *if and only if* the basis is `(π A)^{-1/2} · P̃_n(µ)` without the internal µ-weight; but since the basis of §III.F.1 carries `(Ω·n)` in its inner product, the reflection drops to `diag(1)` in the "total current" basis and the (2n+1) migrates into the basis prefactor. **Test this before changing `reflection_marshak`.** The safe incremental change is step 2 only — fix the µ-weight in the primitives, keep the reflection operator diagonal-Gelbard, and check.

## Verification tests to add after fix

- **Reciprocity (167c)**: compute W^{ρν}_{outer,inner} and W^{νρ}_{inner,outer}, verify `A_outer · W^{ρν}_{oi} = A_inner · W^{νρ}_{io}` with transposed mode indices. Must hold to machine precision.
- **Rank-1 reduction**: set n_modes = 1, verify bit-exact match with F.4 scalar rank-1 Mark closure (already works with Model A/B; must still work after fix).
- **σ_t scan**: the σ_t·R-dependent residual in `peierls_rank_n_closure_not_basis.md` must **collapse to a constant** (ideally zero, or a slow N-convergent residual) after the µ-weight fix.
