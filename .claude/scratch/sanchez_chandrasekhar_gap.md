# Sanchez/Chandrasekhar Green's Functions — Implementation Gap Analysis

**Author:** literature-researcher sub-agent
**Date:** 2026-05-03
**Scope:** Identify literature gaps for implementing Meanings (1) and (2) in ORPHEUS

---

## Headline finding

**The local literature is sufficient for both Meaning (1) and Meaning (2).
No new acquisitions needed.** Furthermore, ORPHEUS already implements
~70% of Meaning (1) for slab + sphere inside `fn_method/.../flux_reconstruction.py`
(KLL 1974 Fredholm-iteration recipe, `_X_negative_real_axis`,
`A(ν)`/`B(ν)` continuum coefficients, `slab_angular_flux_from_scalar` and
`sphere_angular_flux_from_scalar`). What's missing is **packaging discipline**
(promotion to a top-level callable Green's function `G(τ,τ';μ,μ')`) and
**Meaning (2)'s closed-form spectral kernel evaluator** (Sanchez 1986
Eq. (A6) sphere / Pomraning-Siewert 1982 Eq. (21) sphere) as a separate
module that bypasses ray tracing entirely.

This is the OPPOSITE of an acquisition problem — it's a refactoring +
exposure problem. The math is already in the code; the user-facing API
just doesn't surface it as Meanings (1) and (2) yet.

---

## Q1 — Chandrasekhar (Meaning 1) gap

### What Meaning (1) requires constructively

For a 1D plane-parallel slab on `[-a, a]` with isotropic scattering and
vacuum BC, the canonical Chandrasekhar angular Green's function
`G(τ,τ';μ,μ')` is built by:

1. **Discrete eigenpair** `(±ν₀, φ₀±(μ))`: dispersion root + Case discrete
   eigenfunctions.
2. **Continuum eigenfunctions** `φ_ν(μ) = (cν/2)·P[1/(ν-μ)] + λ(ν)δ(μ-ν)`
   for `ν ∈ (-1,1)`.
3. **Wiener-Hopf X-function** for finite-medium half-range bi-orthogonality
   (the Inönü 1973 / Case-Zweifel 1967 form).
4. **Half-range expansion coefficient** `A(ν)` solved by Fredholm iteration
   (KLL 1974 Eqs. 6-7).
5. **Forward angular flux** `ψ(τ, μ) = a₀φ₀(μ)e^{-τ/ν₀} + ∫₀¹ A(ν)φ_ν(μ)e^{-τ/ν}dν`
   for arbitrary internal source.

### Gap analysis (slab + sphere)

| Component | Status | Module | Gap |
|---|---|---|---|
| Dispersion root `ν₀` | DONE | `fn_method/core/dispersion.py::case_nu0`, `case_method/core/dispersion.py::case_atalay_u0` | Both isotropic and anisotropic (Atalay) variants. No gap. |
| Case discrete `φ₀±(μ)` | IMPLICIT | Used inside `slab_angular_flux_from_scalar` via inverse Laplace of slab integral form. | Not exposed as a callable `phi_0_plus(mu, c, nu0)` primitive. |
| Case continuum `φ_ν(μ)` | IMPLICIT | Same — used inside characteristic integration. | Not exposed as a callable. |
| Wiener-Hopf X-function (isotropic) | DONE | `fn_method/slab/flux_reconstruction.py::_X_negative_real_axis` | Private (underscore-prefixed). Needs promotion to public callable. |
| Wiener-Hopf X-function (anisotropic) | DONE | `case_method/core/x_function.py::atalay_X_function` | Atalay Eq. 40 closed-form, scipy + mpmath backends. PUBLIC. Anisotropic c>1 only. No isotropic-c<1 fallback path. |
| McCormick-Kuščer bi-orthogonality | DONE-VIA-ATALAY | `case_method/origins/derivations.py` (per Atalay Eqs. 28-31). | Not exposed as a stand-alone proof; lives inside the criticality solver. |
| Mika 1961 completeness | TEXTUAL | Cited in code/docs, not implemented as a separate object. | Completeness is a theorem, not a numeric primitive — gap is in **documentation**, not code. |
| `A(ν)` Fredholm iteration (slab) | DONE | `fn_method/slab/flux_reconstruction.py::solve_kll_slab_continuum_coefficient` | KLL 1974 Eqs. 6-7. PUBLIC. |
| `B(ν)` Fredholm iteration (sphere) | DONE | `fn_method/sphere/flux_reconstruction.py::solve_kll_sphere_continuum_coefficient` | KLL 1974 Eqs. 13-15. PUBLIC. |
| `ψ(z, μ)` slab from `A(ν)` | DONE | `slab_angular_flux_from_scalar` | Computes the angular flux at any `(z, μ)` via characteristic integration. |
| `ψ(r, μ)` sphere from `B(ν)` | DONE | `sphere_angular_flux_from_scalar` | Same architecture. |
| **Top-level `G(τ,τ';μ,μ')` callable** | **MISSING** | — | The piece that makes it look like "a Chandrasekhar angular Green's function" rather than "a coefficient that you then multiply by something". A 5-arg call: `(c, τ, τ', μ, μ')` returns the response of the medium to a δ-source at `(τ', μ')`. |
| **Sphere → slab parity reduction** (Mitsis Ch.IV) | DONE-IMPLICITLY | `case_method/sphere/one_group.py` uses Atalay's antisymmetric-slab parity reduction. | Mitsis 1963 Ch.IV.5 spherical transformations are not exposed as a generic "sphere Green's function = (1/r)·slab Green's function with parity flip" callable. |
| Cylinder | PARTIAL | `singular_eigenfunction/cylinder/one_group.py` (Westfall-Metcalf 1973) | Returns critical radius and `Φ'(μ)` boundary moment. Does NOT reconstruct `ψ(r, μ)` interior. WM-72 has the formulas (their Eqs. 30-31); ORPHEUS hasn't built them. |

### Required references for Meaning (1) — all LOCAL

| Reference | Local status | Provides |
|---|---|---|
| Mitsis 1963 ANL-6787 | LOCAL — full PhD thesis 160 pp. | Foundational singular-eigenfunction expansion in 3D (Ch.II), plane (Ch.III), sphere (Ch.IV including angular distribution at boundary Fig.7 + at origin Fig.6), cylinder (Ch.V). Mitsis Eq. 4.2-2 = sphere reflection theorem. |
| Atalay 1997 ProgNuclEng 31, 229-252 | LOCAL — full paper | Linear-anisotropic generalization with reflective BC. Eqs. 26+40 for anisotropic X-function. Eqs. 38, 51 for K_j, L_j moments. Eqs. 47-54 sphere = slab antisymmetric parity. Already implemented end-to-end. |
| Kaper-Lindeman-Leaf 1974 NSE 54, 94 | LOCAL | Fredholm iteration for `A(ν)` slab + `B(ν)` sphere; Eqs. 6-7 (slab) and 13-15 (sphere). Critical thicknesses Tables. Already implemented. |
| Westfall-Metcalf 1973 NSE 52, 1-11 | LOCAL | Cylinder singular-eigenfunction expansion. Eqs. 30-32 are the criticality + boundary-moment equations. ORPHEUS has them. Interior reconstruction needs the analogous follow-on. |
| Burkart-Ishiguro-Siewert 1976 | LOCAL | Two-region slab with anisotropic scattering — relevant for multi-region extensions of Meaning (1). |
| Dahl-Sjöstrand 1979 | LOCAL | Eigenvalue spectrum of multiplying slabs and spheres with anisotropic scattering. Useful for verification; not strictly required. |
| Grandjean-Siewert 1979 NSE 69-2 | LOCAL | F_N method numerical results (used by `fn_method`). |
| Siewert-Benoist 1979 NSE 69-1 | LOCAL | F_N method theory. |
| Valougeorgis 1985 PhD | LOCAL | F_N method skeleton + BGK extension (referenced for completeness). |

**No external acquisition required for Meaning (1).**

### What `case_method/sphere/one_group.py` does NOT do (vs. what Meaning (1) requires)

It computes **the criticality determinant** — the value of c such that
the homogeneous problem (zero source) has a non-trivial solution. It
does NOT:

- Evaluate `ψ(r, μ)` for arbitrary internal source.
- Return the point-source response `G(r, r'; μ, μ')`.
- Expose the X-function or the half-range moments as user-callable objects.

These three capabilities ARE present in `fn_method/sphere/flux_reconstruction.py`
(KLL 1974 path — `solve_kll_sphere_continuum_coefficient` +
`sphere_angular_flux_from_scalar`), but they're isotropic-scattering only
and structured around the bare-critical problem (zero source, find R_c)
rather than the fixed-source Green's function.

The gap is therefore **bridge code** — refactor the KLL slab + sphere
paths so the same machinery accepts an arbitrary `q(τ)` source and emits
`(ψ(τ, μ), φ(τ))` as outputs. Math is identical; just exposes the
already-existing `solve_kll_*_continuum_coefficient` with a non-zero
homogeneous-equation RHS.

---

## Q2 — Sanchez 1986 (Meaning 2) gap

### What Meaning (2) requires constructively

For a 1D sphere with linearly anisotropic scattering, specular + diffuse
BC, the source-reduction kernel `K(r, r')` is built by:

1. **Spectral integration over μ** of the within-medium angular Green's
   function (already done by Sanchez 1986 in closed form — Appendix
   Eqs. (A1)-(A7)).
2. **Modified kernel construction** `t(r' → r₀) = t̄(r' → r₀) + t_h(r' → r₀)`
   absorbing the BC reflectivity coefficient `α` (specular) and `β` (diffuse).
3. **Sphere kernel reduction**: `ḡ_2(ρ' → ρ) = (ρ'/ρ)[ḡ_0(ρ' → ρ) - ḡ_0(-ρ' → ρ)]`
   from Eq. (5) — the sphere kernel is the **antisymmetrized slab kernel**.
4. **Slab building block**: `ḡ_0(ρ' → ρ) = (1/2)·E_1(Σ|ρ - ρ'|)` from
   Eq. (6) — the textbook E_1 Bickley-Naylor brick.
5. **Anisotropic flux kernel** `h_α(ρ' → ρ)` for `Σ_{s1}·J(ρ')` term —
   uses `∂_ρ' ḡ_α` plus parity (Eq. (7)).
6. **BC-absorbed kernels** `g̃` and `h̃` in Eq. (13) — closed form via
   Eq. (A6) of the appendix:
   - `g_h(ρ'→ρ) = 2α·∫_{μ₀}¹ T(μ_-)·μ̃⁻¹·cosh(ρμ)·cosh(ρ'μ̃)·exp(-2aμ_-)dμ + (β/(1-βχ_⋆))·Q(ρ')·P(ρ)`
   - where `T(μ) = 1/(1 - α·exp(-τ(μ)))` is the multi-bounce closure,
     `χ⋆ = ∫₀¹ T(μ)χ(μ)μ exp(-2aμ)dμ`, and `P, Q, R` are
     1D-quadrature pieces over `T(μ)·χ(μ)·cosh(ρμ)·exp(-aμ)` etc.
7. **Final flux equation** Eq. (13): linear Fredholm equation for `φ_R(ρ)`.

### Gap analysis (sphere)

| Component | Status | Module | Gap |
|---|---|---|---|
| `E_1(τ)` slab building block | DONE | `derivations/common/kernels.py::E1` | Standard `scipy.special.exp1`. No gap. |
| Sphere kernel = antisymmetrized slab kernel (Eq. 5) | DONE-VIA-VARIANT-α | `peierls_greens_function/greens_function.py` ray-traced. | Variant α uses ray tracing, NOT Eq. (5)'s closed-form parity reduction. The closed-form spectral evaluator is missing. |
| `T(μ) = 1/(1 - α·exp(-τ(μ)))` multi-bounce closure | DONE | `peierls_greens_function/variant_alpha_core.py` | Implemented via geometric series ray tracing. The closed-form `T(μ)` evaluator exists but is **embedded** inside Variant α — could be lifted to `derivations/common/`. |
| Closed-form `g_h(ρ'→ρ)` Eq. (A6) | **MISSING** | — | Variant α evaluates this kernel by ray tracing (per `peierls_greens_function/origins/specular/`) but the **direct quadrature evaluator** of Eq. (A6) doesn't exist. This is the headline gap for Meaning (2). |
| Closed-form `h_h(ρ'→ρ)` Eq. (A6) (anisotropic) | **MISSING** | — | Same. The `2α·∫T(μ_-)·cosh(ρμ)·sinh(ρ'μ̃)·exp(-2aμ_-)dμ` form. |
| `χ⋆`, `K⋆`, `L(x)` auxiliaries (PS-1982 Eqs. 17, 18, 20) | **MISSING** | — | Single 1D quadratures over `T(μ)·χ(μ)·μ·exp(-2Rμ)`. Trivial to add. |
| `F_1(r,x)`, `F_2(r,x)` BC-correction kernels (PS-1982 Eqs. 23-25) | **MISSING** | — | These are the **specular + diffuse** corrections to the bare `E_1(|r-x|) - E_1(r+x)` Sanchez Eq. 5 sphere kernel. They are 1D quadratures over `T(μ)`, `cosh`, `π(x,r,μ) = sqrt(x² - r²(1-μ²))`. |
| Final Fredholm operator on `φ_R(ρ)` Eq. (21) PS-1982 / Eq. (13) Sanchez 1986 | NOT-AS-MEANING-(2) | Variant α does this via ray-traced operator. | Direct closed-form spectral assembly is the deliverable. |
| Linear-anisotropy moment recursion | NEEDED | — | The `J_R(ρ')`, `Q_R(ρ')`, `P_R(ρ)` building blocks of Eq. (13). |
| Cylinder | NOT-IN-PAPER | — | Sanchez 1986 Table 1 lists `α=1` for cylinder but the paper does NOT cover cylinder. PS-1982 is sphere-only. The cylinder analog needs Carlvik 1965 Geneva paper (LOCAL — user-acquired 2026-05-03) and Knyazev 1993 *Atomic Energy* 74(5):368-374, DOI [10.1007/BF00844623](https://doi.org/10.1007/BF00844623) (LOCAL — user-acquired 2026-05-03). **Cylinder Meaning (2) acquisition gap is now CLOSED.** [Earlier drafts cited a phantom "Knyazev-Selivanov 2014" / "Knyazev 2014 JCP" — that reference does not exist; CrossRef shows only the Knyazev 1993 sole-author paper. Phantom hallucinated in `peierls_greens_extensions_lit.md` and propagated; corrected 2026-05-03.] |

### Required references for Meaning (2) — all LOCAL for sphere

| Reference | Local status | Provides |
|---|---|---|
| Sanchez 1986 TTSP 15, 333-343 | LOCAL — read in full | Sphere with linearly anisotropic scattering + specular/diffuse BC. Eq. (1a) generic 1D form. Eq. (5) sphere=antisymmetric slab. Eq. (6) slab E_1 brick. Eq. (13) final Fredholm. **Eqs. (A6)-(A7) = the deliverable**: closed-form `g_h`, `h_h`, `G_h` in 1D quadrature over `T(μ)`, `χ(μ)`, `cosh`, `sinh`. |
| Pomraning-Siewert 1982 JQSRT 28, 503-506 | LOCAL — read in full | Isotropic precursor. Eq. (21) is the desired sphere integral equation `r·I(r) = r·G(r) + ∫₀^R x[(ω/2)I(x) + Q(x)][E_1(|r-x|) - E_1(r+x) + α·F_1(r,x) + β·F_2(r,x)]dx`. Eqs. (22)-(25) give G(r), F_1, F_2 in 1D quadratures. **Direct implementation roadmap.** |
| Mitsis 1963 ANL-6787 | LOCAL | Reflection theorem (Sanchez Eq. 5 derivation), ANL classical sphere = antisymmetric slab. Section IV.5 spherical transformations. |
| Erdmann-Siewert 1968 J.Math.Phys. 9, 81 | NOT-LOCAL | Cited by PS-1982 Ref. [4] as a source for "the integral formulation can be solved in terms of a pseudo-slab problem." | Acquisition wish (low priority — PS-1982 is self-contained). |
| Wu-Siewert 1975 ZAMP 26, 637 | NOT-LOCAL | Cited by PS-1982 Ref. [5] for similar pseudo-slab solution route. | Acquisition wish (low priority). |
| Case-deHoffmann-Placzek 1953 (Introduction to the Theory of Neutron Diffusion, Vol.1, USGov) | NOT-LOCAL (book) | Cited by PS-1982 Ref. [6] as the source of the "integration back along propagation ray" technique used in their Eq. (4). Standard textbook trick — present in many later sources. | Acquisition wish (low priority). |

**Bottom line for Meaning (2) sphere: Sanchez 1986 + PS-1982 + Mitsis 1963 are 100% sufficient.** The gap is implementation, not literature.

### Key implementation facts (extracted from Sanchez 1986 + PS-1982)

- **Slab kernel building block** (PS-1982 doesn't quite use this; Sanchez 1986 Eq. 6 makes it explicit):
  `ḡ_0(ρ'→ρ) = (1/2)·E_1(Σ|ρ-ρ'|)` for the homogeneous-medium isotropic-scattering slab.
- **Sphere parity reduction** (Sanchez 1986 Eq. 5):
  `ḡ_2(ρ'→ρ) = (ρ'/ρ)·[ḡ_0(ρ'→ρ) - ḡ_0(-ρ'→ρ)]`
  i.e. sphere kernel = (ρ'/ρ) × antisymmetric slab kernel. This is the
  **closed-form spectral kernel** for vacuum BC.
- **Specular + diffuse extension** (PS-1982 Eq. 21 / Sanchez 1986 Eq. 13):
  The bare kernel `E_1(|r-x|) - E_1(r+x)` is augmented by `α·F_1(r,x) + β·F_2(r,x)`
  where `F_1, F_2` involve `T(μ)·cosh(...)·exp(-2Rμ_0(...))` 1D quadratures.
- **Multi-bounce closure** `T(μ) = (Re^{d/μ} - e^{-d/μ})/(Re^{-d/μ} - e^{d/μ})`
  is the **same object** that ORPHEUS already evaluates inside Variant α
  via ray tracing — but in Meaning (2) it's a closed-form scalar function
  of μ, not a series.
- **Linear-anisotropic correction term** `Σ_{s1}·J_R(ρ')` is added to
  the source via Sanchez Eq. (13)'s `h̃(ρ'→ρ)·Σ_{s1}·J_R(ρ')` integrand.
  Same `T(μ)` machinery, additional `sinh(ρ'μ̃)` factor and `μ̃⁻¹`
  vs `cosh(ρμ̃)` distinction (Sanchez Eq. A6).

---

## Q3 — Verification matrix sketch

### Three structural paths to `φ(r)` for the bare-critical sphere

| Path | What it computes | How | Where in ORPHEUS |
|---|---|---|---|
| **(α) Variant α / Peierls Green's function** | `φ(r)` via ray-traced multi-bounce angular Green's function | First-leg + bounce-period 2-quadratures + geometric `T(μ)` series | `peierls_greens_function/greens_function.py` |
| **(β) Sanchez kernel** | `φ(r)` via direct closed-form spectral kernel | `E_1` + `F_1`/`F_2` + `χ⋆`/`K⋆`/`L(x)` 1D quadratures | **MISSING — to build** |
| **(γ) Chandrasekhar angular** | `ψ(r, μ)` via singular eigenfunctions, then `∫dμ` | Dispersion root + KLL Fredholm `A(ν)` + `_X_negative_real_axis` | `fn_method/sphere/flux_reconstruction.py::sphere_scalar_flux_kll` |

### Structural-independence boundary

Per `vv-principles` L11: "two derivations are STRUCTURALLY independent
only if they exercise a different integrand or a different identity."

| Pair | Structurally independent? | Reasoning |
|---|---|---|
| (α) vs (β) | **Partially** — same `E_1` building block, different assembly | Both use `E_1`/Bickley-Naylor as the 1D radial brick. (α) assembles via ray-traced multi-bounce; (β) via spectral integration `∫dμ T(μ) cosh(...)`. The **identity** that connects them — that the ray-traced geometric series equals the spectral integral — is provable, not coincidental, so a numerical match is genuine evidence but doesn't probe the `E_1` evaluator itself. |
| (α) vs (γ) | **Strongly independent** | (α) integrand is `q(r')·exp(-Σs)/r²` along characteristics. (γ) integrand is `A(ν)·B_ν(c, R)·sinh(R/ν)/(R/ν)` over Case ν spectrum. Different domain (`r'` vs `ν`), different identity (Peierls vs singular-eigenfunction completeness). A match is the strongest verification claim available. |
| (β) vs (γ) | **Strongly independent** | (β) integrand is `T(μ)·E_1(...)/μ` over surface `μ`. (γ) integrand is `A(ν)·...` over `ν`. Different domain, different identity. The Mika-McCormick-Kuščer completeness theorem connects them, but the numerical paths exercise distinct quadratures. |

### Verification matrix

```
              | (α) Variant α | (β) Sanchez kernel | (γ) Chandrasekhar
--------------|---------------|---------------------|-------------------
verifies (α)  |       —       |  WEAK (shared E_1)  |  STRONG
verifies (β)  | WEAK (shared) |          —          |  STRONG
verifies (γ)  |    STRONG     |       STRONG        |       —
```

**Key insight**: (γ) Chandrasekhar is the **anchor** — it's structurally
independent of both (α) and (β). The (α)↔(β) cross-check is weaker
because both use `E_1` as the 1D building block.

For an L1 verification claim under `vv-principles`, the canonical
combination is:

1. (γ) Chandrasekhar gives the reference number (anchor).
2. (α) and (β) match each other (consistency between two different ways
   of evaluating the same integral kernel).
3. Both (α) and (β) match (γ) (structural independence achieved via
   the singular-eigenfunction route).

This is **stronger than a Sood 2003 benchmark** because the latter is
a single number; here we get three numerical paths converging on the
same answer through distinct integrands.

### Where this matters most

For `c → 1⁺` (deep multiplying limit), (γ) struggles because `ν₀ → ∞`
and the singular-eigenfunction expansion becomes ill-conditioned. (α)
and (β) remain stable because they don't materialize `ν₀`. Conversely,
for very thick multipliers (R ≫ 1 mfp), (α)'s ray-traced bounce series
loses digits to cancellation; (γ) is the cleanest because it's
exponential-in-`R/ν₀`. (β) is intermediate.

So the three paths COVER each other's weak regions — making the
verification claim robust across parameter space.

---

## Q4 — Literature inventory

### Already local (sufficient for both Meaning (1) and Meaning (2) sphere)

| Citation | Local file | What it provides |
|---|---|---|
| Sanchez 1986 TTSP 15, 333-343 | `Sanchez(1986)Integral form...sphere with linearly anisotropic scattering.pdf` | **Meaning (2) headline reference.** Eqs. (1a)-(13) main body + Eqs. (A1)-(A7) appendix. Sphere with linear anisotropy + specular + diffuse BC. |
| Pomraning-Siewert 1982 JQSRT 28, 503-506 | `Pomraning Siewert (1982) On the integral form of the equation of transfer for a homogeneous sphere.pdf` | **Meaning (2) precursor.** Eq. (21) is the explicit isotropic Fredholm form `r·I(r) = r·G(r) + ∫₀^R x[(ω/2)I(x) + Q(x)][E_1(|r-x|) - E_1(r+x) + α·F_1(r,x) + β·F_2(r,x)]dx`. Eqs. (22)-(25) give G, F_1, F_2 in closed-form 1D quadratures over `T(μ)`. |
| Mitsis 1963 ANL-6787 | `Mitsis(1963) Transport Solutions to the Monoenergetic Critical Problems.pdf` | **Meaning (1) foundational.** PhD thesis 160 pp. Ch.II eigenfunctions, Ch.III plane (singular integral for ψ(x,μ), Eq. 3.2 + ψ(x,μ) reconstruction Eq. 3.3d), Ch.IV sphere (sphere = antisymmetric slab), Ch.V cylinder (basis for WM-72). Includes Fig. 6 angular distribution at origin and Fig. 7 angular distribution at boundary — the **point-source response Mitsis already plotted in 1963**. |
| Atalay 1997 ProgNuclEng 31, 229-252 | `Atalay(1997)the reflected slab and sphere criticality problem with anisotropic scattering...pdf` | Meaning (1) anisotropic slab + sphere criticality. Atalay Eq. 40 X-function (Wiener-Hopf form). Eqs. 28-31 parallel-bi-orthogonality (extending McCormick-Kuščer). Eqs. 47-54 sphere = antisymmetric slab. **All implemented in `case_method/`.** |
| Kaper-Lindeman-Leaf 1974 NSE 54, 94 | `Kaper-Lindeman-Leaf(1974)Benchmark Values for the Slab and Sphere...pdf` | KLL 1974 Fredholm-iteration recipe. Eqs. 6-7 slab `A(ν)`, Eqs. 13-15 sphere `B(ν)`. **Already implemented in `fn_method/{slab,sphere}/flux_reconstruction.py`.** |
| Westfall-Metcalf 1973 NSE 52, 1-11 | `Westfall-Metcalf(1972)Singular Eigenfunction Solution...for Finite Radially Reflected Critical Cylinders.pdf` | Cylinder Meaning (1) machinery. Eqs. 30-32. Already implemented for criticality in `singular_eigenfunction/cylinder/one_group.py`; interior `ψ(r,μ)` reconstruction TBD. |
| Burkart-Ishiguro-Siewert 1976 | `Burkart-Ishiguro-Siewert(1976)Neutron transport in two dissimilar media with anisotropic scattering.pdf` | Two-region slab + anisotropic. Useful if Meaning (1) is extended to multi-region. |
| Dahl-Sjöstrand 1979 | `Dahl-Sjostrand(1979)eigenvalue spectrum of multiplying slabs and sphered for monoenergetic neutrons with anisotropic scattering.pdf` | Eigenvalue spectrum diagnostic — verification target for Meaning (1). |
| Grandjean-Siewert 1979 NSE 69-2 | `Grandjean-Siewert(1979)The FN Method in Neutron-Transport Theory. Part II_ Applications and Numerical Results.pdf` | F_N method numerical results. |
| Siewert-Benoist 1979 NSE 69-1 | `Siewert-Benoist(1979)The FN Method in Neutron-Transport Theory. Part I_ Theory and Applications.pdf` | F_N method theory. |
| Garcia 2017 sphere PN | `Garcia(2017) A PN particular solution for the radiative transfer equation in spherical geometry.pdf` | Stable P_N route — a **fourth structurally-independent path** if needed. |
| Garcia 2019/2021 sphere PN | `Garcia(2019)*.pdf`, `Garcia(2021)Accurate spherical harmonics solutions...multi-region spherical geometry.pdf` | Multi-region sphere P_N for fixed-source — useful for extending verification matrix into multi-region. |
| Sood-Foster-Parsons 2003 LA-13511 | `Sood Foster Parsons (1999)Analytical Benchmark Test Set for Criticality Code Verification.pdf`, `Sood-Foster-Parsons(2003)*.pdf` | 75-problem analytical benchmark suite — anchor truth set for both Meanings. |
| Valougeorgis 1985 PhD | `Valougeorgis(1985)The Fn method in kinetic theory PhD thesis.pdf` | F_N method skeleton (BGK extension) — useful methodology cross-reference. |
| Hébert 2009 Ch.3 | `Hebert(2009)Chapter3.pdf` | Modern textbook treatment of Peierls integral form, kernels, and BC absorption — useful context for Meaning (2). |
| Stamm'ler 1983 Ch.4 + Ch.6 | `Stammler(1983)Chapter4.pdf`, `Stammler(1983)Chapter6.pdf` | CP method reference + interface currents. Less load-bearing for Meanings (1) and (2) but useful background. |

### Need from user — STRICTLY OPTIONAL

The local corpus is sufficient for both Meanings on slab + sphere.
The optional acquisitions below are **nice-to-have** for completeness
and citation depth — none are blockers.

| Citation | Why it would help | How to acquire | Priority |
|---|---|---|---|
| Chandrasekhar 1950 *Radiative Transfer* (Dover reprint) | Original H-function formulation, infinite-medium half-space. Useful for documentation cross-reference and as the textbook citation for the H-function nonlinear integral equation. | Likely on user's bookshelf or via library. Out of copyright (Dover reprint exists). | LOW — Mitsis 1963 already cites all needed Chandrasekhar formulas. |
| Case-Zweifel 1967 *Linear Transport Theory* | Foundational textbook for singular eigenfunction theory. Cited everywhere. | Out-of-print but widely available used. Some chapters scanned online. | LOW — Mitsis + Atalay together cover everything we need. |
| Inönü 1973 NSE | Inönü's paper on H/X-functions for finite media. Sanchez 1986 doesn't cite it; Atalay does, indirectly via McCormick-Kuščer 1965. | NSE archive. | VERY LOW — derivation already in Atalay 1997. |
| Mika 1961 NSE | Original completeness theorem for Case eigenfunctions in finite media. | NSE archive. | LOW — referenced in Atalay 1997. The theorem is mathematical, not a numeric primitive — reading the paper is for narrative completeness, not implementation. |
| Case 1960 Annals of Physics 9, 1-23 | Original Case singular-eigenfunction paper. Foundational. | Open-access via Elsevier? Or via citations. | LOW — Mitsis 1963 cites and reproduces all needed material. |
| McCormick-Kuščer 1965 | Half-range bi-orthogonality with linear anisotropy. The **technical foundation** of Atalay's Eqs. 18-21. | NSE archive. | LOW — Atalay 1997 reproduces the relations needed. |
| Kuščer 1955 | Original H-function variant with non-conservative scattering. | Older Slovenian or German journal — harder to find. | VERY LOW — historical context only. |
| Erdmann-Siewert 1968 J.Math.Phys. 9, 81 | Pseudo-slab solution route for sphere. Cited by PS-1982 Ref. 4. | J.Math.Phys. archive — paywall, but JMP allows author preprint. | LOW. |
| Wu-Siewert 1975 ZAMP 26, 637 | Same pseudo-slab route. Cited by PS-1982 Ref. 5. | ZAMP archive — paywall. | LOW. |
| Carlvik 1965 Geneva Conf. Vol.2 p.225 | The **cylinder Wigner-Seitz CP foundational paper**. Cited by Sanchez 1977 Ref. 27. ORPHEUS has Carlvik 1967 (cuboids/finite cylinders) but NOT the 1965 Wigner-Seitz paper. | UN/IAEA archive scan; Argonne report mirror? Worth checking. | MEDIUM — only relevant if the user wants to extend Meaning (2) to **cylinder** geometry (Sanchez 1986 sphere-only). |
| ~~Knyazev-Selivanov 2014 JCP~~ → **Knyazev 1993 (corrected)** | Closed-form 1-D quadrature for the homogeneous linearly-anisotropic cylinder kernel in Bickley-Naylor of order 2+k — cylinder analog of Sanchez 1986 sphere Eq A6. **LOCAL (user-acquired 2026-05-03)**. DOI [10.1007/BF00844623](https://doi.org/10.1007/BF00844623). | The "2014 JCP" citation was a hallucinated phantom; the real paper is *Atomic Energy* 74(5):368-374, May 1993, sole author A. P. Knyazev. |

### Summary verdict: Sub-deliverable A applies

**Local literature is sufficient for both Meaning (1) and Meaning (2)
on slab + sphere with isotropic OR linearly anisotropic scattering.**
No external acquisitions are required to ship working modules.

For **cylinder extension of Meaning (2)** (which Sanchez 1986 does NOT
cover), the user would need to acquire Carlvik 1965 and/or
Knyazev 1993 (corrected from phantom "Knyazev-Selivanov 2014" — that reference does not exist; the real paper is *Atomic Energy* 74(5):368-374, DOI 10.1007/BF00844623, sole author Knyazev, 1993). This is a separate scoping question, not a
blocker for the slab/sphere work.

---

## Implementation order recommendation

### Phase A: Promote existing machinery (1-2 sessions, no new math)

1. **Refactor `fn_method/slab/flux_reconstruction.py` to expose a
   fixed-source path.** Currently the KLL solver assumes the bare-critical
   problem (zero source, finds R_c). Add a parameter `q(τ)` and let the
   continuum coefficient `A(ν)` fit a non-zero RHS.
   - **Deliverable**: `solve_slab_meaning1(c, a, q_func)` returning
     `(phi(z), psi(z, mu))`.
   - **Verification**: with `q_func ≡ 0` and `c = c_critical(a)`,
     should reproduce existing critical-thickness behavior.

2. **Same refactor on sphere**: `solve_sphere_meaning1(c, R, q_func)`
   returning `(phi(r), psi(r, mu))`.

3. **Promote `_X_negative_real_axis` and the discrete/continuum
   eigenfunctions to a public `case_eigenfunctions.py` module** —
   exposes `phi_0_plus(mu, c, nu0)`, `phi_continuum(mu, nu, c)`,
   `X_isotropic(mu, c, u0)` as user-callable primitives.

4. **Build `G_chandrasekhar(c, a, tau, tau_prime, mu, mu_prime)`** as
   the public Meaning (1) entry point. Internally uses (1)-(3).
   - One callable. 5-6 args. Returns one number per call. The
     Chandrasekhar angular Green's function as the user expects.

### Phase B: Build the Sanchez kernel (Meaning 2) for sphere (1-2 sessions)

5. **Implement PS-1982 Eq. (21)** as a separate module
   `peierls_sanchez_kernel/sphere.py`:
   - `g_sphere_isotropic(r, x, R, c)` — the bare kernel
     `(c/2)·x·[E_1(|r-x|) - E_1(r+x)]` (Sanchez Eq. 5/6 + (ω/2) factor).
   - `F1(r, x, R, alpha)` — PS-1982 Eq. 23-24 specular correction.
   - `F2(r, x, R, beta, chi_func)` — PS-1982 Eq. 25 diffuse correction.
   - `solve_sphere_meaning2(c, R, alpha, beta, q_func, chi_func)` —
     the assembled Fredholm operator solved by Nyström quadrature
     in `r ∈ [0, R]`. Returns `phi(r)`.

6. **Extend to Sanchez 1986 linearly anisotropic** (`f_1 ≠ 0`):
   - Add `g_h(r', r)`, `h_h(r', r)` evaluators per Eq. (A6).
   - Add `J_R(ρ')` reduced-current solver per Eq. (9)-(10).
   - Final form: Eq. (13) Fredholm.

### Phase C: Three-way verification matrix (1 session)

7. **Build the cross-check harness**:
   - Generate fixed-source sphere benchmarks (e.g., point source at
     `r' = R/2` with `µ' = 0.5`, evaluate response at `r = R/4`).
   - Run all three paths: (α) Variant α, (β) Sanchez kernel, (γ)
     Chandrasekhar.
   - Tabulate residuals. Should agree to ~1e-8 absolute.

8. **Sphinx page** (`docs/theory/greens_functions.rst`):
   - Document the three Meanings.
   - Show the verification matrix.
   - Cite Mitsis 1963, Sanchez 1986, PS-1982, Atalay 1997, KLL 1974
     with equation numbers.

### Phase D (optional): Cylinder extension

9. **Extend Meaning (1) cylinder** by completing the WM-72 interior
   reconstruction (Eqs. 30-31 with non-zero source). No new literature.

10. **Extend Meaning (2) cylinder** — this is the only path that
    requires Carlvik 1965 + Knyazev 1993 [DOI 10.1007/BF00844623, sole author A. P. Knyazev — earlier drafts misnamed this as "Knyazev 2014" / "Knyazev-Selivanov 2014" which does not exist; corrected 2026-05-03]. Ask the
    user before starting.

### What unblocks what

```
Phase A (Meaning 1 slab + sphere)
    │
    ├── Phase B (Meaning 2 sphere) ─┐
    │                                │
    │                                ├── Phase C (3-way verification)
    │                                │
    └────────────────────────────────┘

Phase D (cylinder) is independent and optional.
```

Phase A is the prerequisite for Phase C because the verification
matrix needs (γ) as the anchor. Phase B is independent of Phase A
mathematically but should follow because Phase A builds the
mathematical scaffolding (X-function, Case eigenfunctions) that the
PS-1982 / Sanchez 1986 evaluator needs at integration nodes.

---

## Notation cross-reference

| ORPHEUS variable | Sanchez 1986 | PS-1982 | Atalay 1997 | Mitsis 1963 | KLL 1974 |
|---|---|---|---|---|---|
| sphere radius R | `a` | `R` | `d` (sphere=R) | `R` | `R` |
| optical depth `τ = Σ_t r` | `ρ` | `r` (optical) | `x` (slab) / `R` (sphere) | `τ` | `r` |
| direction cosine | `μ` | `μ` | `μ` | `μ` | `μ` |
| secondaries-per-collision | `ω_0` (=Σ_s/Σ_t) | `ω` (= Σ_s/Σ_t) | `c` | `c` | `c` |
| linear-anisotropy moment | `ω_1` | (isotropic only) | `f_1` (mean cosine) | (isotropic only) | (isotropic only) |
| specular reflectivity | `α` | `α` | `R` (refl. coeff) | (vacuum only) | (vacuum only) |
| diffuse reflectivity | `β` | `β` | (specular only) | (vacuum only) | (vacuum only) |
| within-medium multibounce | `T(μ)` Eq. (A4) | `T(μ)` Eq. (14) | `T(R,μ)` Eq. 33, `T_1(R,μ)` Eq. 50 | (none) | `T(R,μ)` |
| Wiener-Hopf X-function | (not used) | (not used) | `X(μ)` Eq. 26/40 | `X(z)` (Sec. 3.2-2B) | `X(-ν)` (= our `_X_negative_real_axis`) |
| dispersion root | `ν₀` (implicit) | `μ`-spectrum continuous | `ν₀ = i u₀` | `ν₀` (real for c<1, imaginary for c>1) | `ν₀ = i u₀` |
| sphere kernel | `g̅_2(ρ'→ρ)` Eq. (5) | `[E_1(|r-x|) - E_1(r+x)]` from Eq. (21) | (uses parity reduction) | spherical kernel Sec. 4.4-1 | (used implicitly) |
| BC-correction kernels | `g_h, h_h` Eq. (A6) | `F_1, F_2` Eqs. 23-25 | `T(R,μ)`, `K_j` moments | (vacuum) | (vacuum) |

The notation is consistent enough that **no notation conflicts** block
implementation. The biggest gotcha is that Sanchez 1986 uses `α=0,1,2`
as a geometry index (slab/cyl/sphere) while everywhere else `α` is the
specular reflectivity. ORPHEUS code already uses `alpha` for
reflectivity, so the geometry index would need to be renamed (e.g.,
`geom_index` or `dim`).

---

## Caveats and notation gotchas to flag during implementation

1. **`α` overload**: Sanchez 1986 Table 1 uses `α ∈ {0, 1, 2}` for
   geometry. Sanchez 1986 main body uses `α ∈ [0, 1]` for specular
   reflectivity. Pick one, stick with it. ORPHEUS convention: `alpha`
   = specular reflectivity ∈ [0,1].

2. **`ω` vs `c`**: Sanchez 1986 + PS-1982 use `ω = Σ_s/Σ_t`. Atalay
   1997 + KLL 1974 + ORPHEUS use `c = ν̄·Σ_f / Σ_t + Σ_s/Σ_t` — i.e.
   `c` includes fission. For pure-scattering tests `c = ω`. **For
   multiplying problems** the literature `c` is the relevant
   "secondaries-per-collision" parameter (which is what ORPHEUS uses).

3. **`τ(μ)` chord-length convention**: Sanchez 1986 Eq. (A2) uses
   `τ(μ) = 2aμ` for the SPHERE (where `a` is the radius). PS-1982
   Eq. (14) writes `T(μ) = 1/[1 - α·exp(-2Rμ)]`. Same identity:
   the chord crossing a homogeneous sphere along axial direction
   `μ` has length `2Rμ`. **For an off-axis chord** the length is
   `2R√(1 - (h/R)²)` where `h` is the impact parameter — this is
   the `π(x, r, μ)` function in PS-1982 Eq. (10) and Sanchez Eq. (A5).

4. **Sphere reflection theorem sign**: Mitsis 1963 / Sanchez 1986 Eq. 5:
   `ḡ_2(ρ'→ρ) = (ρ'/ρ)·[ḡ_0(ρ'→ρ) - ḡ_0(-ρ'→ρ)]`. The minus sign is
   load-bearing. Atalay 1997 Eq. 47 writes `ψ(x,μ) = -ψ(-x,-μ)` (sphere
   = odd modes). The relationship: on the antisymmetric slab `[-R, R]`,
   the sphere flux corresponds to the **odd-mode** subspace.

5. **`χ⋆` saturation**: PS-1982 Eq. (17) and Sanchez 1986 Appendix
   define `χ⋆ = ∫₀¹ T(μ)·χ(μ)·μ·exp(-2Rμ)dμ`. The factor
   `(1 - β·χ⋆)⁻¹` in Eq. (16)/Eq. (13) can be SINGULAR if `β·χ⋆ → 1`
   — that's the diffuse-BC-only super-critical limit. For `β = 0`
   (pure specular) it never matters. For `β > 0` and `χ ≡ 1/π`
   (Lambertian) it's a numerical-stability flag.

6. **`A(ν)` Fredholm-iteration starting guess**: KLL 1974 starts the
   iteration from the F_N approximation. ORPHEUS already does this
   (per `solve_kll_slab_continuum_coefficient`). Don't re-derive.

7. **Atalay X-function is c>1 only**: `case_method/core/x_function.py`
   raises `ValueError` for `c <= 1`. The KLL X-function
   `_X_negative_real_axis` works for any `c`. **For Meaning (1) the
   KLL form is the right one to promote** — Atalay's adds anisotropy
   but is restricted to multiplying media.

---

## Summary table

| Question | Answer |
|---|---|
| Q1 (Meaning 1 gap) | NO new acquisitions needed. Code gap = expose existing KLL machinery as a fixed-source `G(τ,τ';μ,μ')` callable; promote `_X_negative_real_axis` to public; build `case_eigenfunctions.py` exposing `φ₀±` and `φ_ν`. |
| Q2 (Meaning 2 gap) | NO new acquisitions needed for sphere. Sanchez 1986 + PS-1982 + Mitsis 1963 are sufficient. Code gap = build `peierls_sanchez_kernel/sphere.py` evaluating PS-1982 Eq. (21) directly via 1D quadratures. |
| Q3 (Verification matrix) | (γ) Chandrasekhar is the structurally-independent anchor. (α) and (β) share `E_1` so their cross-check is weaker. The TRIPLE check (α)≈(β)≈(γ) is L1-grade evidence per `vv-principles`. |
| Q4 (Literature) | Sub-deliverable **A** — already-local literature is sufficient for both Meanings on slab + sphere. Cylinder Meaning (2) requires Carlvik 1965 (LOCAL 2026-05-03) + Knyazev 1993 *Atomic Energy* 74(5):368-374, DOI 10.1007/BF00844623 (LOCAL 2026-05-03). [Earlier drafts cited a phantom "Knyazev 2014" that does not exist; corrected.] |

---

## Pointers for the main agent

If the user accepts this plan, the natural sub-agent dispatch sequence is:

1. **method-implementer** for Phase A (refactor + expose Meaning 1).
   Brief should include: `fn_method/{slab,sphere}/flux_reconstruction.py`
   path, KLL 1974 references, the L1 acceptance criterion (zero-source
   reduces to existing criticality answer).

2. **method-implementer** for Phase B (build Sanchez kernel sphere).
   Brief should include: PS-1982 Eq. (21)-(25) as the equation of
   record, Sanchez 1986 Eq. (A6) as the anisotropic extension,
   `derivations/common/kernels.py::E1` as the building block.

3. **test-architect** for Phase C (verification matrix). Brief should
   include: the three paths' expected agreement at 1e-8 absolute,
   regions where each path may degrade (high-c, large-R, near-vacuum).

4. **archivist** for Phase C documentation. Brief should include: the
   three meanings, the local PDFs for each, the verification matrix.

The **critical first action** before any implementation is the user's
explicit confirmation that:

(a) The naming of `peierls_sanchez_kernel/` (or alternative name like
    `sanchez_kernel/`) is acceptable architecturally.
(b) Phase D (cylinder) is deferred unless they want Carlvik 1965
    acquired now.
(c) The `α` overload issue is resolved one way (recommendation:
    `alpha` = specular reflectivity, geometry index = `geom`).
