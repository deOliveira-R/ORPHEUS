# TODO: Fix Cylindrical 1D SN Azimuthal Diamond-Difference

## Status: BROKEN for heterogeneous problems

The cylindrical 1D SN sweep produces **exact eigenvalues for homogeneous problems** (all groups, both quadratures) but **wildly wrong and divergent eigenvalues for heterogeneous problems** (keff oscillates with mesh refinement instead of converging).

Spherical 1D SN works correctly for both homogeneous and heterogeneous.
Cartesian SN works correctly for both.

## Architecture

### Files involved

| File | Role |
|------|------|
| `02.Discrete.Ordinates/sn_sweep.py` | `_sweep_1d_cylindrical()` — the broken function |
| `02.Discrete.Ordinates/sn_geometry.py` | `SNMesh._setup_cylindrical()` — computes per-level α |
| `02.Discrete.Ordinates/sn_quadrature.py` | `ProductQuadrature`, `LevelSymmetricSN` — provide level structure |
| `tests/test_sn_cylindrical.py` | Tests — 20 pass, 1 xfail (heterogeneous) |

### How the cylindrical sweep works

The 1D cylindrical transport equation (per unit axial height) is:

```
(η/r) ∂(rψ)/∂r  −  (1/r) ∂(ξψ)/∂φ  +  Σ_t ψ  =  Q / sum(w)
```

where:
- η = sin(θ)cos(φ) — radial direction cosine (= `mu_x`)
- ξ = sin(θ)sin(φ) — azimuthal direction cosine (= `mu_y`)
- μ = cos(θ) — axial direction cosine (= `mu_z`)
- A = 2πr — cylindrical face area per unit height
- V = π(r²_out − r²_in) — cylindrical cell volume per unit height

The angular quadrature is organized by **μ-level**: on each level (fixed μ_z), there are M azimuthal ordinates at angles φ_0, φ_1, ..., φ_{M-1}.

**Azimuthal redistribution**: the `−∂(ξψ)/∂φ` term couples successive azimuthal ordinates on each level. The standard discretization uses:

```
α_{m+1/2} = α_{m-1/2} + w_m · ξ_m
```

with boundary conditions α_{1/2} = 0, α_{M+1/2} = 0 (by azimuthal symmetry of ξ).

**The discrete balance equation** (integrated over cell i, ordinate m):

```
η_m [A_{i+1/2} ψ_{i+1/2,m} − A_{i-1/2} ψ_{i-1/2,m}]
  − [α_{m+1/2} ψ_{i,m+1/2} − α_{m-1/2} ψ_{i,m-1/2}]
  + Σ_t V_i ψ_{i,m}
  = Q V_i / sum(w)
```

**NOTE THE MINUS SIGN** before the α terms. This is the key difference from spherical, which has a plus sign.

### Current implementation (in `_sweep_1d_cylindrical`)

The sweep currently uses `|α|` (absolute values) with a **positive** sign convention — identical to the spherical sweep. This is **wrong** for the cylindrical equation which has a negative sign.

The DD closures are:
- Spatial: `ψ_{i+1/2} = 2ψ − ψ_{i-1/2}` (standard diamond-difference)
- Angular: `ψ_{m+1/2} = 2ψ − ψ_{m-1/2}` (standard diamond-difference)

After substitution, the current code solves:
```
ψ = [QV/sum_w + |η|(A_in+A_out)·ψ_in + (|α_out|+|α_in|)·ψ_angle] / [2|η|A_out + 2|α_out| + Σ_t V]
```

**This has the WRONG sign on the α terms** — it should be minus, not plus.

## The problem

### Symptom

| Test | Result |
|------|--------|
| Homogeneous 1G/2G/4G, Product quad | **Exact** to machine precision |
| Homogeneous 1G/2G/4G, LS S4 quad | **Exact** to machine precision |
| Heterogeneous 1G, 5 cells/region | keff = 1.14 (CP ref: 0.99) |
| Heterogeneous 1G, 10 cells/region | keff = 0.88 (CP ref: 1.00) |
| Heterogeneous 1G, 20 cells/region | keff = 0.47 (CP ref: 1.00) |
| Heterogeneous 1G, 40 cells/region | keff = 0.22 (diverging!) |
| Heterogeneous 2G, Product 4×8 | keff = 0.52 (CP ref: 0.74) |
| Heterogeneous 2G, Product 8×8 | keff = 0.90 (oscillates!) |

Homogeneous works because the spatial flux is flat — no net streaming, so the redistribution term vanishes (αψ cancels in the balance). The bug only manifests when there are spatial gradients.

### Root cause: sign convention

The cylindrical azimuthal redistribution has a **minus** sign:

```
−[α_{m+1/2} ψ_{m+1/2} − α_{m-1/2} ψ_{m-1/2}] / V
```

After DD substitution `ψ_{m+1/2} = 2ψ − ψ_{m-1/2}`:

```
−[α_{m+1/2}(2ψ − ψ_{m-1/2}) − α_{m-1/2} ψ_{m-1/2}] / V
= −2α_{m+1/2}ψ/V + (α_{m+1/2} + α_{m-1/2})ψ_{m-1/2}/V
```

The correct DD equation is:
```
denom = 2|η|A_out − 2α_out + Σ_t V
numer = QV/sum_w + |η|(A_in+A_out)ψ_in − (α_out+α_in)·ψ_angle
```

Note: **signed α values** (not |α|), with minus sign.

### Why the correct sign doesn't work

When α_{m+1/2} is positive and large (middle of the azimuthal sweep), the denominator `2|η|A_out − 2α_out + Σ_t V` can become **zero or negative**, causing overflow. This happens particularly at inner cells where A_out is small and α is large.

The α values for a Product(2,4) quadrature are:
```
Level 0: α = [0.00, 0.00, 1.28, 1.28, 0.00]
```

At α_out = 1.28: `denom = 2(0.82)(0) − 2(1.28) + 0.5(0.79) = −2.16`, which is negative.

### What the spherical sweep does differently

In the spherical sweep (`_sweep_1d_spherical`), the angular redistribution has a **positive** sign:

```
+[α_{n+1/2} ψ_{n+1/2} − α_{n-1/2} ψ_{n-1/2}] / V
```

After DD substitution:
```
denom = 2|μ|A_out + 2|α_out| + Σ_t V
```

The denominator is always positive (all terms positive), so the scheme is unconditionally stable. This is why spherical works but cylindrical doesn't — the opposite sign makes the denominator potentially negative.

## Solutions attempted

### 1. Use |α| with positive sign (current code)
**Result**: Stable but wrong for heterogeneous. Exact for homogeneous.
**Problem**: Wrong sign convention; effectively treats cylindrical like spherical.

### 2. Use signed α with correct minus sign
**Result**: Overflow/NaN — denominator goes negative.
**Problem**: Standard DD is not positive-definite with the correct sign.

### 3. Skip angular DD update at level boundaries (α_out = 0)
**Result**: No improvement — same wrong keff.
**Problem**: The issue is at interior ordinates (large α), not boundaries.

### 4. Clamp angular face flux to non-negative (np.maximum)
**Result**: No improvement — same keff values as pure DD.
**Problem**: The face fluxes aren't negative; the ψ values themselves are wrong.

### 5. Step method for angular variable (ψ_{m+1/2} = ψ)
**Result**: Converges with refinement (1.21→1.13→1.05) but slowly (first-order in angle). Still ~5% error at 20 cells/region. 2G gives 0.49 vs CP 0.74.
**Problem**: First-order accuracy in angle is insufficient.

### 6. Zero redistribution (α = 0 manually)
**Result**: Smooth spatial convergence (1.38→1.39→1.39). Wrong keff (no angular physics) but stable.
**Problem**: Proves the spatial streaming is correct; the redistribution is what breaks things.

## Hypothesis for the fix

The standard solution in the nuclear engineering literature is the **Lewis & Miller starting-direction treatment** (§4.5.4 in "Computational Methods of Neutron Transport", 1984). The key ideas:

### Track αψ instead of ψ

Instead of updating `psi_angle = 2ψ − psi_angle` (which is ψ_{m+1/2}), track the **product** `Φ = α·ψ`:

```
Φ_{m+1/2} = α_{m+1/2} · ψ_{m+1/2}
```

The balance equation in terms of Φ:
```
η streaming − [Φ_{m+1/2} − Φ_{m-1/2}]/V + Σ_t V ψ = QV/sum_w
```

With the DD closure `ψ = ½(ψ_{m-1/2} + ψ_{m+1/2})`:

```
ψ_{m+1/2} = Φ_{m+1/2} / α_{m+1/2}  (when α ≠ 0)
ψ_{m-1/2} = Φ_{m-1/2} / α_{m-1/2}  (when α ≠ 0)
ψ = ½(Φ_{m-1/2}/α_{m-1/2} + Φ_{m+1/2}/α_{m+1/2})
```

At the boundaries (α = 0), Φ = 0 regardless of ψ, so the boundary condition is naturally satisfied.

### The "starting direction" pseudo-ordinate

At the first ordinate on each level (α_{1/2} = 0), the incoming angular flux is determined by conservation rather than DD extrapolation. Lewis & Miller defines a "starting direction" flux:

```
ψ_{start} = (total source into this direction) / (total removal from this direction)
```

This provides the initial Φ_{1/2} for the azimuthal sweep.

### Reference

- Lewis, E.E. and Miller, W.F. Jr., *Computational Methods of Neutron Transport*, §4.5 "Curvilinear Coordinates", Wiley, 1984.
- Morel, J.E. and Montry, G.R., "Analysis and Elimination of the Discrete Ordinates Flux Dip", *Transport Theory and Statistical Physics*, 13:5, 1984.
- Alcouffe, R.E. et al., "PARTISN: A Time-Dependent, Parallel Neutral Particle Transport Code System", LA-UR-05-3925, Los Alamos, 2005.

## Testing regime for the fix

Run these in order. Each must pass before proceeding to the next:

### 1. Unit test: αψ product conservation
Verify that Σ (α_{m+1/2} ψ_{m+1/2} − α_{m-1/2} ψ_{m-1/2}) over all m on a level = 0 (telescoping sum with α boundary = 0).

### 2. Single-cell uniform source
Two cells, 1G, Σ_t=1, Q=1. After many sweeps, volume-averaged φ → Q/Σ_t = 1.

### 3. Homogeneous 1G exact
20 cells, S4 product quadrature. Must match analytical k_inf to < 1e-6.

### 4. Homogeneous 2G/4G exact
Same mesh. Must match analytical k_inf.

### 5. Heterogeneous 1G convergence
Fuel (r<0.5) + moderator (r<1.0), 5/10/20/40 cells per region.
keff must converge monotonically with refinement.
Compare with CP reference (should differ by ~1-5% from white-BC approximation).

### 6. Heterogeneous 2G resolution independence
Product(4×8) and Product(8×8) must agree to < 5%.
Product(4×8) and Product(8×16) must agree to < 3%.

### 7. Particle balance heterogeneous
Production / absorption = keff for 2G fuel+moderator.

### 8. Cross-check with CP cylindrical
1G heterogeneous: SN and CP should agree to ~5% (white-BC vs reflective difference).

### 9. Both quadratures agree
Product and Level-Symmetric S_N (once LS is fixed) must give close keff on heterogeneous.

## Additional known issue: Level-Symmetric S_N structure

The `LevelSymmetricSN` quadrature groups ±μ_z hemispheres on the same level (level 0 has 16 ordinates with both μ_z = +0.41 and μ_z = −0.41). The cylindrical sweep needs each level to have a SINGLE μ_z value. This is a separate bug from the DD sign issue.

**Fix**: In `_build_level_symmetric()`, create separate levels for +μ_z and −μ_z instead of grouping by |μ_z|. Or restructure `level_indices` to split mixed-hemisphere levels.

The `ProductQuadrature` does not have this issue — each level has exactly one μ_z value.

## Files to modify for the fix

1. `02.Discrete.Ordinates/sn_sweep.py` — rewrite `_sweep_1d_cylindrical()` with αψ tracking
2. `02.Discrete.Ordinates/sn_quadrature.py` — fix LS level structure (separate issue)
3. `tests/test_sn_cylindrical.py` — remove xfail once tests pass
