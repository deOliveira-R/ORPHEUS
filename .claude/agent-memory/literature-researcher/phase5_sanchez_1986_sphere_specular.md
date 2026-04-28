---
name: Sanchez 1986 — sphere integral form with specular BC (Phase 5 textbook source)
description: Extracted content from Sanchez 1986 TTSP DOI 10.1080/00411458608210456. Eq. (A6) IS the continuous-µ multi-bounce specular kernel for homogeneous sphere — the textbook root for ORPHEUS Phase 5. T(µ) = 1/(1 − α·exp(−2aµ)) algebraically identical to ORPHEUS expected 1/(1 − exp(−σ·2Rµ)) with α=1, a=ΣR. Multi-region NOT in paper. µ-weight convention vs M1 sketch flagged for SymPy verification.
type: project
---

# Sanchez 1986 TTSP 14 — Sphere integral kernel with specular BC (Phase 5 textbook source)

PDF location: `/workspaces/ORPHEUS/Integral form of the equation of transfer for a homogeneous sphere with linearly anisotropic scattering.pdf`. DOI: 10.1080/00411458608210456. Pages 333-343 (the appendix carries the implementation kernel).

## TL;DR

**Sanchez 1986 IS the textbook continuous-µ multi-bounce specular kernel for the homogeneous sphere.** The "multi-bounce" factor appears as

```
T(µ) = 1 / [1 − α·exp(−τ(µ))]
```

with `τ(µ) = 2aµ` the chord across a homogeneous sphere of optical radius `a` at direction-cosine `µ`, and `α ∈ [0,1]` the **specular reflection coefficient**. This is **algebraically identical** to ORPHEUS's expected `1/(1 − exp(−σ·2Rµ))` once we set `α = 1` (perfect specular) and recognise `a = Σ·R` (Sanchez non-dimensionalises optical length).

The Phase 5 sphere kernel is **literally Eq. (A6) with β = 0, α = 1, ω₁ = 0**.

Two non-trivial gaps for ORPHEUS production:

1. **Multi-region**. Paper is homogeneous-sphere only. Sanchez's closed-form `T(µ)` extends with piecewise τ(µ) under specular BC (every bounce traverses ALL regions, so the chord τ stays the same per bounce); only the source-arc and receiver-arc factors `cosh(ρµ)` need replacement by per-region path attenuation via `optical_depth_along_ray`.
2. **µ-weight convention**. ORPHEUS's M1 sketch has `µ` in the numerator of `f_mb`; Sanchez's `T(µ)` has no `µ` in numerator but the kernel carries `µ_*⁻¹` (in `g_h`) or no µ-weight (in `h_h`). The two forms are equivalent up to the volume-Jacobian; **SymPy verification required before production implementation** (Phase 5a step 1).

## Sanchez's setup

**Geometry**: 1-D plane / cylinder / sphere unified by integer `α` (Sanchez Table 1):

| Geometry | α | a₋ | ρ |
|---|---|---|---|
| Slab | 0 | -a | depth `z` |
| Cylinder | 1 | 0 | polar radius `ρ` |
| Sphere | 2 | 0 | radius `r` |

**Length scale**: total optical distance is `τ`; `Σ` is the absorption coefficient. The paper uses `a` for the *optical* outer radius, and `2aµ` is the optical chord at direction-cosine `µ`. Conversion to ORPHEUS: `a = Σ·R`.

**Phase function** (linearly anisotropic): `f(Ω·Ω') = (1/4π) (ω₀ + ω₁ Ω·Ω')` with `Σ_s = Σ·ω₀`, `Σ_s1 = Σ·ω₁`, `Σ_a = (1 − ω₀)·Σ`.

**Boundary condition** (Eq. A3.a):

```
I_in(µ) = α·I_out(µ) + β·χ(µ)·∫₀¹ I_out(µ') µ' dµ'
```

- `α ∈ [0,1]` is the specular reflection coefficient.
- `β ∈ [0,1]` is the diffuse reflection coefficient.
- `χ(µ)` is the diffuse re-emission profile, normalised `∫₀¹ χ(µ) µ dµ = 1`.
- `α + β = 0` recovers vacuum BC.
- `α = 1, β = 0` is **perfect specular reflection** (the ORPHEUS Phase 5 case).
- `α = 0, β = 1` is isotropic diffuse reflection (Lambertian with `χ = const`).
- `α + β = 1` is conservative (no leakage).

## Key equations (verbatim)

### (1a) — Unified 1-D integral transport equation

```
I(ρ) = ∫_{a₋}^{a} [g_α(ρ' → ρ)·F(ρ') + h_α(ρ' → ρ)·Σ_s1·J(ρ')] dρ' + G(ρ)
```

ORPHEUS map: `I → φ`, `Q → Q`, `g_α → K_vol_vol`, `h_α → 0` (isotropic), `G → K_bc`. **`G(ρ)` IS the Phase 5 kernel.**

### (3) — Green's function defining BVP

```
(Ω·∇ + Σ) t = δ(r − r') δ_2(Ω·Ω'),    a₋ ≤ ρ < a,
t(r',Ω' → r,Ω) = α·t(r', Ω' → r, Ω_R) + (β/2π)·χ(µ) ∫_{2π+} t(r',Ω' → r,Ω") · Ω"·n dΩ",   |ρ| = a, Ω·n ≤ 0
```

with `Ω_R = Ω + 2µ·n` the specularly reflected direction. The BC line is Sanchez's "homogeneous BC incorporated into the kernel" — the modified kernel that the abstract refers to.

### (A4) — Closed-form bounce sum (THE PHASE 5 KERNEL ROOT)

```
I_in(µ) = (1 / 2πA) · e^{−τ_+} · T(µ_+) · [α · δ(µ − µ_+) / µ_+ + β/(1 − β·χ_*) · χ(µ) · T(µ)]
```

with

```
T(µ) = 1 / [1 − α·e^{−τ(µ)}]                          (THE multi-bounce factor)
χ_*  = ∫₀¹ T(µ) · e^{−τ(µ)} · χ(µ) · µ dµ              (diffuse self-coupling)
```

For homogeneous sphere of optical radius `a`:

```
τ(µ) = 2aµ          (chord across the sphere at cosine µ)
τ_+  = aµ_+ − ρ'µ'  (interior point ρ' → surface)
τ_-  = aµ_- + ρµ    (surface → interior point ρ)
T(µ) = 1 / [1 − α·e^{−2aµ}]
```

**This is IDENTICALLY the Phase 5 kernel** with `α = 1`.

### (A5) — Reflected-streaming Green's function

```
t_h(r' → r) = (1 / 2πA) · e^{−τ_+ − τ_-} · T(µ_+)
              · [α · δ(µ_- − µ_+) / µ_+ + β/(1 − β·χ_*) · T(µ_-) · χ(µ_-)]
```

Specular leg `α · δ(µ_- − µ_+) / µ_+` preserves direction cosine; diffuse leg `β/(1−β·χ_*) · T(µ_-) · χ(µ_-)` re-emits with profile `χ`. Both legs multiplied by the "trapped-photon" factor `T(µ_+)·exp(−τ_+ −τ_-)`.

### (A6) — Multi-bounce volume-volume kernels (THE PHASE 5 IMPLEMENTATION)

For the homogeneous sphere:

```
g_h(ρ' → ρ) = 2α ∫_{µ₀}^{1} T(µ_-) · µ_*⁻¹ · cosh(ρµ) · cosh(ρ'µ_*) · e^{−2aµ_-} dµ
              + (β / (1 − β·χ_*)) · Q(ρ') · P(ρ)

h_h(ρ' → ρ) = 2α ∫_{µ₀}^{1} T(µ_-) · cosh(ρµ) · sinh(ρ'µ_*) · e^{−2aµ_-} dµ
              + (β / (1 − β·χ_*)) · R(ρ') · P(ρ)
```

with auxiliary integrals:

```
P(ρ)  = ∫₀¹ T(µ_-) · χ(µ_-) · cosh(ρµ) · e^{−aµ_-} dµ
Q(ρ') = (2 ρ'² / a²) ∫₀¹ T(µ_+) · cosh(ρ'µ') · e^{−aµ_+} dµ'
R(ρ') = (2 ρ'² / a²) ∫₀¹ T(µ_+) · sinh(ρ'µ') · e^{−aµ_+} µ' dµ'
```

and

```
µ₀² = max(0, 1 − (ρ'/ρ)²)                 (chord-visibility cone)
µ_*² = ρ'² − ρ²(1 − µ²)                    (chord-projection identity)
```

**This is the ORPHEUS Phase 5 sphere kernel verbatim.** Setting `ω₁ = 0` drops `h_h`; the term in `T(µ_-)` is the multi-bounce factor; the `cosh(ρµ)·cosh(ρ'µ_*)` structure is the spherical analog of the slab `cosh` Green's function (it comes from the even-extension of `J_R(ρ)` through `ρ = 0`).

The `e^{−2aµ_-}` factor is the chord attenuation across the *full* sphere; `T(µ_-)` then sums the geometric series over the bounces.

### (A7) — Boundary-source kernel `G(ρ)` (incident `I_0` ≠ 0)

```
G(ρ) = 4π [ ∫₀¹ T(µ_-) · I_0(µ_-) · cosh(ρµ) · e^{−aµ_-} dµ
            + (β / (1 − β·χ_*)) · K* · P(ρ) ]
```

with `K* = ∫₀¹ T(µ') · I_0(µ') · e^{−2aµ'} · µ' dµ'`.

This is the kernel ORPHEUS multi-region extension reuses: in a multi-region sphere, the inner region sees the outer region's outgoing surface flux as its incident `I_0(µ)`.

## Continuous-µ kernel form vs ORPHEUS M1 sketch

ORPHEUS expected:

```
K_bc^mb,sph(r_i, r_j) = 2 ∫₀¹ G_in(r_i, µ) · F_out(r_j, µ) · [µ / (1 − e^{−σ·2Rµ})] dµ
```

Sanchez's Eq. (A6) **is exactly this**, with the explicit identification:

| ORPHEUS variable | Sanchez variable | Notes |
|---|---|---|
| `G_in(r_i, µ)` | `cosh(ρ·µ) e^{−aµ_-}` | "in-streaming" factor at receiver |
| `F_out(r_j, µ)` | `cosh(ρ'·µ_*) µ_*⁻¹` (g) or `sinh(ρ'·µ_*)` (h) | source-side outgoing factor |
| `µ / (1 − e^{−σ·2Rµ})` | `T(µ_-) = 1/(1 − α e^{−2aµ_-})` | multi-bounce factor (`α=1`) |
| `σ·2R` | `2a` | optical diameter (Sanchez non-dimensionalises) |
| Outer factor `2` | `2α` | leading `α` appears at the kernel level |

**µ-weight subtlety**: Sanchez's `T(µ)` does NOT carry a `µ` prefactor inside; the surrounding kernel has `µ_*⁻¹` (in `g_h`) or no µ-weight (in `h_h`). ORPHEUS M1 sketch's `µ/(1−exp(−σ·2Rµ))` is `T(µ)·µ` — equivalent under the volume-Jacobian but the prefactor placement is different. **SymPy verification needed before implementing** (Phase 5a step 1).

The integration limit `µ₀ = √max(0, 1 − (ρ'/ρ)²)` is the **chord-visibility cone**: at receiver radius `ρ` and source radius `ρ'`, the trajectory must clear the inner sphere of radius `min(ρ,ρ')`, so `µ ≥ µ₀`. ORPHEUS plumbing for this exists already in `geometry.rho_max`.

## Multi-region extension

**Sanchez does NOT treat multi-region.** The closed-form `T(µ)` in Eq. (A4) holds because every bounce sees the same chord `2a`. In a multi-region sphere with regions `R_1, …, R_N` and absorption `Σ_1, …, Σ_N`, the bounce sees a piecewise τ(µ) = Σ_k Σ_k·ℓ_k(µ).

**Critical observation for specular BC**: under perfect specular reflection (α=1), every bounce traverses ALL regions in the same order with the same chord lengths (the geometry is identical bounce-to-bounce), so the multi-bounce factor `T(µ) = 1/(1 − e^{−τ(µ)})` extends DIRECTLY with piecewise τ(µ). The closed-form bounce sum HOLDS.

What changes is the source-arc and receiver-arc factors `cosh(ρµ)` and `cosh(ρ'µ_*)`, which carry the optical depth `σ·ρµ` of the partial chord; in multi-region these become piecewise integrals via `optical_depth_along_ray`.

Two extension paths flagged by the agent (one is more conservative than the analysis above):

**(i) Outermost-shell-only treatment.** Treat the outer shell with Sanchez's closed form, propagate the surface-flux through the outer-to-inner interface as `I_0(µ)` in Eq. (A7).

**(ii) Direct multi-region piecewise τ(µ)** (per critical observation above). Replace `e^{−2aµ}` with `e^{−τ(µ)}` and `cosh(ρµ)` etc. with the multi-region path attenuation. Closed form survives.

ORPHEUS Phase 5 plan adopts path (ii). Path (i) is fallback.

## Anisotropic scattering treatment

Sanchez handles linearly anisotropic via `h_h` kernel mapping net current `J(ρ)` (not energy `I(ρ)`) into emission. For ORPHEUS today (`ω₁ = 0`), `h` drops entirely. When ORPHEUS adds linearly anisotropic scattering later, port both `g_h` and `h_h` from Eq. (A6); structure is parallel.

## Numerical recommendations

Sanchez **does not state a numerical method** in this paper — it is a theoretical reduction. The integrals are 1-D proper integrals on `[µ₀, 1]` (specular leg) and `[0, 1]` (diffuse leg).

Integrand behaviour:

- At `µ → 0`, `T(µ) → 1/(1 − α)` finite for `α < 1`. For `α = 1` (perfect specular), `T(µ) ~ 1/(2aµ)` — an integrable simple pole. Sanchez notes p. 339 that the `µ₀ = √max(0, 1 − (ρ'/ρ)²)` cutoff resolves this when `ρ' < ρ`.
- Smooth on the open interval — Gauss-Legendre converges spectrally except near the µ → 0 pole at α = 1.
- For `α → 1` the µ → 0 simple pole requires either a change of variables (`µ = u²` or `µ = sin t`) or a Gauss-Jacobi rule weighted by `µ`. **THIS is the numerical concern for ORPHEUS Phase 5** — exactly the grazing-µ divergence the Phase 4 matrix-Galerkin form choked on.

Fredholm structure: Eq. (13) is a Fredholm equation of the second kind in `φ_R(ρ)`, smoothing kernel, spectrally compact. Eigenvalue extraction via power iteration on the integral operator works directly — no `(I − T·R)` matrix inversion needed.

## Citation graph

Sanchez 1986 cites:

1. **Pomraning & Siewert (1982)** *J. Quant. Spec. Rad. Transf.* 28, 503 — "Integral form of the equation of transfer for a homogeneous sphere". Isotropic-scattering precursor (`ω₁ = 0` case). PRIORITY: pull as cross-check truth source for ORPHEUS Phase 5 sphere homogeneous.
2. **Chandrasekhar (1950/1960)** *Radiative Transfer* — slab integral-equation methods.
3. **Williams (1971)** *Mathematical Methods in Particle Transport* — slab Fredholm.
4. **Sanchez (1975, CEA N-1831)** — prior CEA report on Galerkin resolution.
5. **Sanchez (1974, CEA N-1793)** — prior CEA Galerkin report; reduction Eq. (12) traces here.
6. **Sanchez & McCormick (1982)** NSE 80, 481 — "A review of neutron transport approximations". Sanchez 1986 uses kernel decomposition Eq. (5) and conservation identity Eq. (6) from there.

Notable absences: NO citation of Davison, Hébert (postdates), Stamm'ler.

## Phase-5-critical pages

- **pp. 333-334**: Setup, BC formulation, motivation.
- **p. 335 (Table 1)**: notation map slab/cyl/sphere unified.
- **p. 337 (Eq. 11, 13)**: reduced integral equation for vacuum and reflective BC.
- **pp. 338-339**: explicit form of the kernels with BC. **Most important page for ORPHEUS implementation.**
- **pp. 340-342 (Appendix, A1-A7)**: derivation of `t_h` and the closed forms `g_h`, `h_h`, `G`. **The implementation roadmap.**

## Relationship to ORPHEUS Phase 5 plan

| Item | Use Sanchez? | Guidance | Gaps |
|---|---|---|---|
| Reduced integral equation | YES — Eq. (1a) + (13) | Direct plug-in; isotropic case sets `ω₁ = 0` and drops `h` | None |
| Multi-bounce factor `1/(1−e^{−σ2Rµ})` | YES — `T(µ) = 1/(1 − α·e^{−2aµ})` in Eq. (A4) | Algebraically identical with `α=1`. Sanchez's non-dim: `a = Σ·R` | µ-weight convention SymPy check |
| Green's function with BC | YES — Eq. (A1), (A5) | Modified kernel `t = t̄ + t_h`. Sanchez's `t_h` IS the Phase 5 BC closure | None for sphere |
| Eigenvalue / k_eff structure | PARTIAL | Power iteration on (13); standard | k-eigenvalue extraction not stated |
| Multi-region | NO — homogeneous-only | Closed form survives via piecewise τ(µ) under specular BC | Implement piecewise τ(µ) plumbing |
| Anisotropic scattering | YES — full | Use `h_h` from Eq. (A6) when `ω₁ ≠ 0` | None — already general |
| Specular vs diffuse parametrization | YES — Eq. (A3.a) | `α + β` with arbitrary `χ(µ)` | None |
| Quadrature method | NO | Gauss-Jacobi or `µ → u²` for the grazing pole at `α=1` | Implementation choice |
| Comparison with Pomraning-Siewert | YES — direct generalization | PS 1982 is `ω₁=0` case. Use as cross-check truth | Pull PS 1982 PDF |

## File index

- PDF location: `/workspaces/ORPHEUS/Integral form of the equation of transfer for a homogeneous sphere with linearly anisotropic scattering.pdf`
- Page count: 11 numbered pages (pp. 333-343) + appendix
- Phase-5-critical pages: pp. 338-342 (kernel + appendix)
