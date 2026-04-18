# Phase A / C empirical investigation — 2026-04-18

## Setup

Bare homogeneous 1G 1-region white-BC eigenvalue problem:
- Σ_t = 1, Σ_s = 0.5, νΣ_f = 0.75, k_∞ = 1.5
- Geometry: {cylinder-1d, sphere-1d}
- Cell radius R ∈ {1, 2, 5, 10} MFP
- Rank N ∈ {1, 2, 3, 5, 8}
- Quadrature: `n_angular = n_rho = n_surf_quad = 24`, `dps = 25`

## Variants tested

**V1 — CURRENT (legacy mode 0 + `(ρ_max/R)²` Jacobian for n ≥ 1 + `R = diag(2n+1)`)**
- Mode 0 routes through legacy `compute_P_esc` / `compute_G_bc` (no Jacobian)
- Modes n ≥ 1 use `compute_P_esc_mode` (with `(ρ_max/R)²` Jacobian)
- Marshak reflection `R = diag(1, 3, 5, …, 2N-1)`

**V2 — CONSISTENT JACOBIAN (`compute_*_mode` for ALL modes + `R = diag(2n+1)`)**
- Mode 0 also uses `compute_P_esc_mode` with Jacobian
- All other structure matches V1

**V4 — CANONICAL COEFFICIENT SPACE (V2 + `R = B^{-1}`)**
- Same integrand as V2
- `R = B^{-1}` where `B_nm = ∫_0^1 μ P̃_n(μ) P̃_m(μ) dμ` (cosine Gram matrix)
- The Marshak closure `J^-_n = J^+_n` (moment matching) is equivalent
  to `α^-_m = α^+_m` (coefficient matching) if `B^{(N)}` is invertible;
  the Gram matrix pivots between moment space and coefficient space.

**V6 — COSINE-WEIGHTED INTEGRAND (`μ_exit · P̃_n` inside + `R = B^{-1}`)**
- `compute_P_esc_mode_v6` uses `∫ W_Ω · μ_exit · P̃_n(μ_exit) · K_esc(τ) dΩ`
- No `(ρ_max/R)²` Jacobian; cosine weight is `μ_exit` directly in integrand
- `R = B^{-1}` (canonical cosine-moment inversion)

## Headline numbers

### Sphere (k_∞ = 1.5)

```
            V1 (current)                V2 consistent             V4 B^{-1}                  V6 cosine-integrand
R    N=1  N=2  N=3  N=5  N=8      N=1  N=2  N=3  N=5  N=8    N=1   N=2   N=3   N=5   N=8     N=1  N=2  N=3  N=5  N=8
1.0  27%  1.2% 2.3% 2.5% 2.5%     50%  29%  25%  24%  24%    29%   16%   15%   15%   15%     1.1% 19%  19%  19%  19%
2.0  7.2% 3.4% 3.8% 3.8% 3.8%     36%  27%  26%  26%  26%    26%   22%   22%   22%   22%     17%  24%  24%  24%  24%
5.0  0.7% 0.3% 0.3% 0.3% 0.3%     15%  15%  14%  14%  14%    14%   14%   14%   14%   14%     11%  11%  11%  11%  12%
10.0 0.3% 0.2% 0.2% 0.2% 0.2%     5.3% 5.2% 5.2% 5.2% 5.2%   5.2%  5.2%  5.2%  5.2%  5.2%    6.3% 6.6% 6.7% 6.8% 6.9%
```

### Cylinder (k_∞ = 1.5)

```
            V1 (current)                V2 consistent             V4 B^{-1}
R    N=1  N=2  N=3  N=5  N=8      N=1  N=2  N=3  N=5  N=8    N=1   N=2   N=3   N=5   N=8
1.0  21%  8.3% 27%  60%  107%     41%  17%  0.5% 28%  65%    23%   2.6%  13%   40%   78%
2.0  7.3% 1.5% 4.9% 10%  17%      27%  20%  16%  11%  6.4%
5.0  2.1% 1.4% 1.2% 0.8% 0.2%     9.9% 9.5% 9.3% 9.0% 8.7%
10.0 1.1% 1.1% 1.0% 1.0% 0.9%
```

## Observations

1. **V1 is empirically the best** overall: sphere converges to 2.5 % plateau,
   cylinder rank-1 matches legacy Mark. Other variants are universally worse
   except V2 cylinder rank-3 at R=1 MFP hits **0.45 %**, the canonical
   DP_2 target — demonstrating that the *structure* of the closure can
   achieve high-rank convergence, but with a different mode-0 convention
   that degrades rank-1 accuracy.

2. **The mode-0 convention is decisive.** V1 (legacy `P_esc` without
   Jacobian) gives the Mark rank-1 answer (20–30 % error depending on R).
   V2/V4/V6 (canonical with Jacobian or cosine weight on mode 0) give
   very different rank-1 answers — sometimes better (V6 sphere R=1: 1 %),
   sometimes worse (V2 sphere R=1: 50 %). There is no single "right"
   mode-0 convention without calibration against a reference.

3. **Changing `R` from `diag(2n+1)` to `B^{-1}`** does not rescue the
   degraded rank-1 convention. V4 (V2 + `B^{-1}`) improves sphere from
   24 % plateau to 15 % plateau — a marginal improvement — but cylinder
   still diverges at high N.

4. **V6's cosine-integrand rank-1 is accidentally good** for sphere at
   R=1 (1 %) but its rank-N never improves. The rank-1 match is not a
   canonical result; it is a geometric coincidence at thin R where
   `μ_exit ≈ 1` for most rays.

## Diagnosis

The "canonical Marshak DP_N" cannot be achieved by a trivial rewrite
of the integrand + `R` matrix within the ORPHEUS V1 convention. Multiple
entangled factors resist simple corrections:

- The mode-0 legacy convention is not the canonical partial-current
  moment; it is the *isotropic-source escape probability*. These agree
  at `r = 0` in a sphere but differ for general `r`.
- The `(ρ_max/R)²` factor in `compute_P_esc_mode` is an
  empirical correction that works for mode n ≥ 1 but does not
  follow from a clean derivation of the cosine-weighted moment.
- The cylinder `compute_G_bc_mode` uses a different angular
  convention (surface-centered Ki_1/d) from `compute_P_esc_mode`
  (observer-centered Ki_2), which compounds the mode-n issues.
- `R = B^{-1}` is the canonical coefficient-space closure but only
  works when both `P` and `G` use matched canonical integrands, which
  requires rewriting both sides consistently.

The canonical answer likely requires calibration against the
[Stepanek1981] slab DP_N k_eff tables to anchor the mode-0
normalization, followed by a consistent rewrite of both `P` and `G`
for all modes. This is more involved than the simple integrand swap
initially hoped; see [Sanchez1982] §III.F.1 for the full DPN setup.

## Status of Phases A and C

- **Phase A** (sphere plateau at 2.5 %): **deferred**. No simple fix
  identified. Likely requires Stepanek calibration.
- **Phase C** (cylinder divergence at N ≥ 3): **deferred**. Same reason.
- **V1** remains the shipped implementation with its known limitations.

## Scripts

- `diag_rank_n_09_phaseAC_scan.py` — V1, V2, V3 variant scan (V3 never
  completed; partial data in V1+V2 output)
- `diag_rank_n_10_V4_canonical_Binv.py` — V4 (canonical `R = B^{-1}`)
- `diag_rank_n_11_V5_cylinder_Ki_kplus2.py` — V5 cylinder Ki_{k+2}
  expansion (did not complete due to `mpmath.quad` for high-k Ki being
  slow; see [Knyazev1993] polynomial expansion for future efficient
  impl)
- `diag_rank_n_12_V6_cosine_integrand.py` — V6 cosine-weighted integrand

## References

- Sanchez, R. & McCormick, N. J. "A Review of Neutron Transport
  Approximations." *Nucl. Sci. Eng.* **80**, 481–535 (1982).
- Stepanek, J. "The DP_N Surface Flux Integral Neutron Transport Method
  for Slab Geometry." *Nucl. Sci. Eng.* **78**, 171–179 (1981).
- Knyazev, A. P. "Solution of the transport equation in integral form in
  a one-dimensional cylindrical geometry with linearly anisotropic
  scattering," *Atomic Energy* **74** (5), 385–389 (1993).
