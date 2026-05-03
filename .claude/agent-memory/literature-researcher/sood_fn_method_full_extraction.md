---
name: Sood/Forster/Parsons LA-13511 (1999) — full benchmark suite extraction
description: Complete map of the 75-problem analytical benchmark test set, naming convention, transport-equation reductions, and the canonical 5-case complexity ramp for F_N reimplementation. Companion to sood_2003_cylinder_benchmarks.md.
type: reference
---

# LA-13511 (1999) — full report extraction

PDF: `/workspaces/ORPHEUS/scratch/literature/Sood Foster Parsons (1999)Analytical Benchmark Test Set for Criticality Code Verification.pdf`
(also at `/tmp/la13511.pdf`, identical MD5 e017ef0c…). 57 pages, 67 tables.
2003 PNE summary paper is the journal-published condensation of this report.

## 1. Report structure overview

The report is **a compendium of published analytical-method results**, NOT a methods paper. Sood et al. tabulate 75 critical configurations whose r_c, k_inf, and selected scalar-flux ratios were computed in peer-reviewed journals using **Case singular eigenfunctions, F_N method, B_N method, and Green's-function method** (p. 2). The report itself contains:

- Sections I-IV: motivation, naming, scope (no math).
- Section V: transport-equation reductions for 1G/2G slab (Eqs. 1-17).
- Sections VI-IX: 67 cross-section + critical-dimension + flux tables, organized by (groups, scattering, geometry, material).
- Appendix A: closed-form k_inf algebra for 1G, 2G, 3G, MG infinite media (Eqs. 18-76).

**Crucial finding**: the F_N method machinery (Case eigenfunctions, Wiener-Hopf factorization, H-functions) is NOT in this report. It lives in the references — primarily:
- Kaper, Lindeman, Leaf (1974) NSE 54, 94 — slab + sphere F_N machinery.
- Westfall, Metcalf (1973) NSE 52, 1 — finite radially-reflected cylinders (singular eigenfunction).
- Grandjean, Siewert (1979) NSE 69, 161 — F_N method canonical Part II.
- Sanchez, Ganapol (1983) NSE 64, 61 — cylinder linearly anisotropic.
- Kornreich, Ganapol (1997) NSE 126, 293 — Green's-function method.
- Sjostrand (1986), Garis (1991), Sahni-Sjostrand papers — reflected-sphere/slab eigenvalue spectra.

To implement the F_N method in SymPy, one of those papers (Kaper-Lindeman-Leaf for slab+sphere; Westfall-Metcalf for cylinder) must be obtained. The LA-13511 report is the **truth set** (problem definitions + reference values), not the **method specification**.

### Naming convention (Table 1, Section III)

```
Bare:       <Material>-<Groups>-<Scattering>-<Geometry>     e.g. PUb-1-0-CY
Reflected:  <Material>-<Reflector(thickness)>-<Groups>-<Scattering>-<Geometry>
                                                        e.g. UD2O-H2O(10)-1-0-CY
```

| Field           | Codes                                                       |
| --------------- | ----------------------------------------------------------- |
| Fissile         | PU, U, UD2O, UAL, URR (+ optional `a`/`b`/`c`/`d` variant)  |
| Reflector       | bare (none), H2O, Fe-Na                                     |
| Reflector(N)    | thickness in **mfp** of the reflector. `IN` = infinite.     |
| Groups          | 1 / 2 / 3 / 6                                               |
| Scattering      | 0 = isotropic, 1 = P_1, 2 = P_2                             |
| Geometry        | IN = infinite, SL = slab, CY = cylinder, SP = sphere, ISLC  |

Total: **75 problems** (43 1G, 30 2G, 1 3G, 1 6G). 24 infinite, 24 slabs, 9 cylinders, 14 spheres, 4 ISLC. The Variant α cylinder benchmark we already use (`Ua-1-0-CY`, problem 13) is in this set with r_c = 1.72500292 mfp / 5.284935 cm.

## 2. Transport-equation reduction chain (LA-13511 Section V — verbatim with ORPHEUS notation)

**Eq. 1** — General k_eff eigenvalue equation:
```
Ω·∇Ψ(r,E,Ω) + Σ_t(r,E)Ψ(r,E,Ω)
   = ∫dE'∫dΩ' Σ_s(r,E'→E,Ω'→Ω) Ψ(r,E',Ω')
   + χ(E) ∫dE' [ν(r,E')/(4π·k_eff)] Σ_f(r,E') ∫dΩ' Ψ(r,E',Ω')
```

**Eq. 2** — 1D slab, azimuthally symmetric (Ω·∇ → µ∂/∂x), continuous in E.

**Eq. 5** — 1G isotropic, slab, in optical thickness z = Σ_t·x:
```
µ ∂Ψ/∂z + Ψ(z,µ) = (c/2) ∫_{-1}^{1} Ψ(z,µ') dµ',     c = (Σ_s + νΣ_f/k_eff)/Σ_t
```

This is the **canonical 1G monoenergetic transport eigenvalue problem** Sood-style. For sphere/cylinder, the streaming term gets the geometry factor (Bell-Glasstone Eq. 1-49). For c > 1 a finite critical r_c exists.

**Eq. 7** — Linearly anisotropic with separable scattering/fission angular kernels:
```
µ ∂Ψ/∂z + Ψ = (c/2) ∫ Π(µ_o, β) Ψ(z,µ') dµ'
Π(µ_o, β) = β·Π_s(µ_o) + (1-β)·Π_f(µ_o),    β = Σ_s / (Σ_s + νΣ_f)
```

**Eq. 10** — 2G isotropic, slab, matrix form:
```
µ ∂Ψ̄/∂z + Σ̄ Ψ̄ = (1/2) C̄ ∫ Ψ̄(z,µ') dµ'
Ψ̄ = [Ψ_1; Ψ_2],  Σ̄ = diag(Σ, 1),  C̄ = [c_11 c_12; c_21 c_22],  c_ij = Σ_ij / Σ_2
```

(Sood uses `2 = fast, 1 = slow` — opposite of ORPHEUS convention.)

**Eq. 15-17** — 2G with linear anisotropy, with c_ij,l = (2l+1)·Σ_ij,l / Σ_2.

### Appendix A — k_inf closed forms (the easy ones — no F_N needed)

- **1G**: k_inf = c·Σ_t·νΣ_f / [(Σ_t-Σ_s)(Σ_s+νΣ_f)]    (Eq. 20)
- **2G**: explicit ratio of Σ-products (Eq. 28); upscatter-free version Eq. 29.
- **3G**: matrix-determinant form (Eq. 42, 59) + worked O'Dell example.
- **MG general**: k = ν̄Σ_f^T (Σ̄_t - Σ̄_s)^{-1} χ̄ ν̄Σ_f φ̄  (Eq. 76, single matrix inversion).

These ALL admit closed-form SymPy verification with no special functions. Branch 1 of the F_N project starts here — they're the trivial cases that exercise infinite-medium algebra.

## 3. The proposed 5-case complexity ramp

Selected for monotone difficulty + structural diversity. Each labeled with `(problem #)` from Tables 2-4.

### Case 1 — `PUa-1-0-IN` (problem 1) — infinite medium, 1G isotropic

| Quantity   | Value                                            |
| ---------- | ------------------------------------------------ |
| Geometry   | infinite, no spatial structure                   |
| XS (Pu-239 a) | ν=3.24, Σ_f=0.0816, Σ_c=0.019584, Σ_s=0.225216, Σ_t=0.32640 |
| c          | 1.50                                             |
| Reference  | k_inf = **2.612903** (constant Ψ everywhere)     |

**Derivation path**: Eq. 20 directly. SymPy: `k = c*Σ_t*ν*Σ_f / ((Σ_t-Σ_s)*(Σ_s+ν*Σ_f))`. Verifiable to machine precision. **MVP starting point.**

### Case 2 — `Ua-1-0-SL` (problem 12) — bare slab, 1G isotropic

| Quantity   | Value                                            |
| ---------- | ------------------------------------------------ |
| Geometry   | bare slab                                        |
| XS (U-235 a) | ν=2.70, Σ_f=0.06528, Σ_c=0.013056, Σ_s=0.248064, Σ_t=0.32640 |
| c          | 1.30                                             |
| Reference  | r_c = **0.93772556 mfp / 2.872934 cm** (Ref. 26 = Kaper-Lindeman-Leaf 1974) |
| Flux       | φ(r/r_c=0.25)=0.9669506, φ(0.5)=0.8686259, φ(0.75)=0.7055218, φ(1.0)=0.4461912 |

**Derivation path**: Kaper-Lindeman-Leaf 1974 NSE 54, 94. The slab F_N is rank-N approximation to the half-range Wiener-Hopf factorization X(z, µ); critical condition is det(matrix-built-from-X-evaluations-at-Case-eigenvalue-points) = 0. **First geometry that requires F_N — the canonical entry point for the method.**

### Case 3 — `Ua-1-0-CY` (problem 13) — bare infinite cylinder, 1G isotropic

| Quantity   | Value                                            |
| ---------- | ------------------------------------------------ |
| Geometry   | infinite cylinder, finite radius                 |
| XS         | same as Case 2 (U-235 a, c=1.30)                 |
| Reference  | r_c = **1.72500292 mfp / 5.284935 cm** (Refs 27,28 = Westfall-Metcalf 1972/1973) |

**This is the cylinder benchmark we already use** — Variant α agrees to 8.5e-6, per `sood_2003_cylinder_benchmarks.md`. The Westfall method uses **singular eigenfunction expansion + Bessel-function angular flux representation** (NOT F_N for cylinder — Westfall uses a structurally distinct technique). Sanchez-Ganapol 1983 NSE 64, 61 (Ref. 31) is the F_N variant for cylinder linearly anisotropic and would also be a Branch 2 target.

### Case 4 — `Ua-1-0-SP` (problem 14) — bare sphere, 1G isotropic

| Quantity   | Value                                            |
| ---------- | ------------------------------------------------ |
| Geometry   | bare sphere                                      |
| XS         | same as Case 2 (U-235 a, c=1.30)                 |
| Reference  | r_c = **2.4248249802 mfp / 7.428998 cm** (Ref. 26 = Kaper-Lindeman-Leaf 1974) |

**Same paper as Case 2** — Kaper-Lindeman-Leaf 1974 derive both slab AND sphere F_N machinery. Sphere is the natural Variant α geometry (Plan 2 already targets sphere) — independent F_N solver here gives the structurally-independent cross-check Plan 2 needs. **High priority.** Note: no flux ratios published for Ua-1-0-SP in this table (Table 14), but `PUb-1-0-SP` does have flux ratios (Table 8: 0.93538006 / 0.75575352 / 0.49884364 / 0.19222603 at r/r_c = 0.25/0.5/0.75/1.0).

### Case 5 — `PU-2-0-IN` (problem 44) — infinite medium, 2G isotropic

| Quantity   | Value                                            |
| ---------- | ------------------------------------------------ |
| Fast group (g=2): | ν_2=3.10, Σ_2f=0.0936, Σ_2c=0.00480, Σ_22s=0.0792, Σ_12s=0.0432, Σ_2=0.2208, χ_2=0.575 |
| Slow group (g=1): | ν_1=2.93, Σ_1f=0.08544, Σ_1c=0.0144, Σ_11s=0.23616, Σ_21s=0.0, Σ_1=0.3360, χ_1=0.425 |
| Reference  | k_inf = **2.683767**, φ_2/φ_1 = **0.675229** |

**Derivation path**: Eq. 28 (no upscatter, χ_2 + χ_1 = 1), plus Eq. 32 for flux ratio. Pure linear algebra — 2x2 matrix determinant + division. **The simplest 2G test of the multi-group machinery.**

**Bonus 6th case** for when MR is added: `PU-2-0-SP` (problem 46) — bare 2G Pu sphere, r_c = 1.15513 mfp / 5.231567 cm (Ref. 8 = Siewert-Thomas 1986 NSE 94, 264). This is the **2G sphere F_N** — the union of Cases 4 + 5 — and the natural target for testing 2G Variant α extensions. Reference 36 (Forster PhD thesis 1970) is the slab equivalent.

## 4. Continuous flux verdict — NO closed-form reconstruction from the report

The report tabulates flux ratios at exactly 4 sample points: `r/r_c ∈ {0.25, 0.5, 0.75, 1.0}` (Tables 8, 14, 18, 23, 42, 55). **These are not closed-form sample points** — they are independently-published numerical evaluations of the Case-eigenfunction expansion at those specific abscissae. The report does not give the eigenfunction coefficients.

**To reconstruct φ(r) at ANY r ∈ [0, R]** requires re-deriving the F_N expansion from the cited paper:
- For slab/sphere: Kaper-Lindeman-Leaf 1974 supplies the X-function recursion + the φ(z) formula in terms of the dispersion-law roots and Case eigenvalue coefficients.
- For cylinder: Westfall-Metcalf 1973 supplies the Bessel-function representation; flux is a finite linear combination of `J_0(z·µ)/(1-µ²)^(1/2)` evaluations integrated over the F_N quadrature. Sanchez-Ganapol 1983 gives the linearly-anisotropic version.

Once you have those expansions in SymPy/numpy, you can evaluate φ(r) on a dense grid and the **4 published points become 4 verification anchors** (5-digit truth values at known r). The flux on the OTHER points is **only as good as the F_N reimplementation** — i.e. you cross-check Variant α against your own F_N curve, with the published 4 points anchoring both.

## 5. Implementation feasibility — where SymPy walls live

The **easy half** — closed-form on inspection, full SymPy:

- **All k_inf cases (Appendix A)**: 1G/2G/3G/MG. Pure rational algebra in cross sections. Branch 1 = SymPy expression manipulation; reference values match to machine precision.
- **Bare-slab eigenvalue problem at finite N** (Kaper-Lindeman-Leaf): the F_N matrix entries involve Case eigenvalues `ν_0, ν_1` (transcendental roots of `1 - c·ν·atanh(1/ν) = 0` for c>1) and integrals of the form `∫ µ^k / (ν₀ ± µ) dµ`. The integrals are **closed-form rational+log** (computable in SymPy). The transcendental root finding has to fall to mpmath `findroot`. Branch 1 = symbolic matrix entries; Branch 2 = high-precision numerical root + det evaluation.
- **2G matrix k_inf**: SymPy `Matrix(...).det()` / `solve` works directly.

The **hard half** — Bessel functions and Wiener-Hopf:

- **Cylinder F_N or singular eigenfunction (Westfall-Metcalf, Sanchez-Ganapol)**: the angular flux involves modified Bessel `K_n` and Bickley-Naylor `Ki_n` integrals. SymPy supports symbolic `besselk` BUT specific integrals like `∫_0^{π/2} cos(...)·K_0(...) dθ` are NOT in SymPy's closed-form table; they fall to numerical quadrature. Branch 1 = symbolic preamble (define the integral, factor the dispersion law), Branch 2 = mpmath quadrature.
- **Wiener-Hopf X-function for half-space**: in slab/sphere F_N, the half-range factorization X(z,µ) requires evaluating a singular Cauchy integral. Closed forms via the Chandrasekhar H-function exist in **integral form** (`X(z,µ) = exp[(1/π)·∫_0^1 atan(...)·(1/(t-µ)) dt]`). SymPy can carry the integral symbolically but cannot evaluate it; use mpmath quadrature with at least 30 dps for 6-digit critical-radius accuracy. The N→∞ asymptote is the truth value; convergence tests at N=8/12/16/24 should give 5+ digits.
- **Multi-region / reflected geometries**: interface conditions become a linear system of size `(N+1)·(num_regions)`. Conditioning degrades for thick reflectors. Use mpmath if scipy.linalg fails to converge.

**SymPy walls** (the user's prior concern):

1. **High-N F_N matrices (N ≥ 20)**: SymPy `Matrix.det()` is symbolic; for 20×20 with rational entries it's tractable but slow. For 40×40 it crawls. Mitigation: keep N small (Sood reports converge at N=8-16) and verify N-convergence numerically with mpmath rather than symbolically.
2. **Special-function integrals (Bickley-Naylor, Bessel-K integrated)**: not in SymPy's table. Mitigation: **two-tier strategy** — SymPy for the structural form (which integrals appear), mpmath for the numerics. This mirrors what Plan 2 Variant α already does.
3. **Transcendental dispersion law roots**: `1 - c·ν·atanh(1/ν) = 0`. SymPy can manipulate the equation symbolically but not solve it. Mitigation: mpmath `findroot` with bisection bracket from c-asymptotics.

**Net feasibility**: Cases 1, 2, 4, 5 fully feasible in (SymPy + mpmath). Case 3 (cylinder) needs Bessel-K numerical quadrature but the rest is symbolic. All five cases can match published values to **at least 5 decimal places**, matching Sood's claimed accuracy.

## 6. Folder structure proposal

Existing convention (per Plan 2 architecture) lives at:
- `orpheus/derivations/peierls_nystrom/`
- `orpheus/derivations/peierls_greens_function/`

These are **method-keyed**, not problem-keyed. Recommend matching:

```
orpheus/derivations/fn_method/
├── __init__.py
├── core/                      # method primitives shared across cases
│   ├── case_eigenvalues.py    # ν_0, ν_1 transcendental roots; X-function
│   ├── moments.py             # ∫µ^k/(ν±µ) integrals (slab moment library)
│   ├── h_function.py          # Chandrasekhar H/X half-range factorization
│   └── bessel_kernels.py      # Ki_n, K_n integrals for cylinder
├── slab/                      # Kaper-Lindeman-Leaf 1974 NSE 54, 94
│   ├── one_group.py           # solve_critical_slab(c) → r_c
│   └── two_group.py           # Forster 1970 / Stewart 1974 thesis
├── sphere/                    # Kaper-Lindeman-Leaf 1974 (same paper as slab)
│   ├── one_group.py
│   └── multi_region.py        # reflected sphere — Sjostrand 1986, Sahni-Sjostrand
├── cylinder/                  # Westfall-Metcalf 1973 NSE 52, 1
│   ├── one_group.py
│   └── linearly_anisotropic.py  # Sanchez-Ganapol 1983 NSE 64, 61
├── multi_group/               # Siewert-Thomas 1986 NSE 94, 264 sphere; Forster slab
│   ├── infinite_medium.py     # Appendix A k_inf cases (trivial)
│   └── two_group_sphere.py    # PU-2-0-SP test
└── benchmarks/
    └── la13511.py             # the 75-problem catalogue, machine-readable
```

Cross-check tests (Branch 2 of `algebra-of-record`) should live alongside the existing peierls tests in `tests/derivations/`, with a new `test_fn_la13511_*.py` set verifying that **Variant α and the F_N reimplementation agree on the same physics to ≥6 digits** for every overlapping case (sphere primary, cylinder bonus, slab if Variant α slab is added).

The Siewert-Thomas 1986 paper (Ref. 8) is cited as the analytic source for **all 2G bare slab + sphere** problems (45-49). Securing that PDF (NSE 94, 264 — likely paywalled, try OpenAlex/HAL/Sci-Hub with DOI 10.13182/NSE86-A18620 or similar) unlocks Cases 5 + bonus.

## 7. Action items for the user

1. **Acquire Kaper-Lindeman-Leaf 1974** (NSE 54, 94 — DOI 10.13182/NSE74-A23308 likely) — unlocks slab + sphere F_N for Cases 2, 4. **Top priority.**
2. **Acquire Westfall-Metcalf 1973** (NSE 52, 1) — already cited in `sood_2003_cylinder_benchmarks.md` as the cylinder "F_N method" source; unlocks Case 3.
3. **Acquire Siewert-Thomas 1986** (NSE 94, 264) — unlocks Case 5 + 2G sphere bonus.
4. **Optional: Sanchez-Ganapol 1983** (NSE 64, 61) — linearly anisotropic cylinder; complements problem 36 (`Ua-1-1-CY`) cross-check.

If any of these are paywalled, the Plan-A literature pull route (HAL, Zenodo, OpenAlex `oa_url`, Semantic Scholar `oa_pdf_url`) should be tried before falling back to alternative formulations. Kaper-Lindeman-Leaf is old enough (1974) that an OSTI lab-report version may exist — Kaper was at Argonne.

## Cross-references

- `sood_2003_cylinder_benchmarks.md` — the 2003 PNE summary; cylinder-only.
- `bickley_naylor_sphere_white_bc.md` — sphere white-BC closure (relevant when extending Case 4 to reflected variants).
- `phase4_cylinder_peierls.md` — cylinder Peierls Ki_n machinery (overlaps with Westfall-Metcalf cylinder F_N).
