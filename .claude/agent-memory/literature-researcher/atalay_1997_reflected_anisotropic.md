---
name: Atalay 1997 — reflected slab+sphere with linearly anisotropic scattering
description: Case singular-eigenfunction method for reflected slab (and sphere via antisymmetric trick) with linearly anisotropic scattering. Primary candidate for unifying the reflected+anisotropic Sood cases.
type: reference
---

**Citation**: M. A. Atalay (1997), *Prog. Nucl. Energy* **31**(3), 229-252.
DOI: `10.1016/0149-1970(95)00094-1`. NOT in Zotero (suggest user add it —
high-value primary source).

# 1. Method classification

**Case singular-eigenfunction method** + **first-order Fredholm
iteration** on the half-range bi-orthogonality machinery of
McCormick–Kušcer 1965. Atalay derives FOUR additional half-range
relations (Eqs. 28–31) parallel to the McCormick–Kušcer set
(Eqs. 18–21) — these are the new technical contribution that lets
the singular-EF method close the reflected-slab problem with
anisotropy (previously stuck on a deficit in the bi-orthogonality
relations; see p.230, last paragraph).

Linearly anisotropic only (Σ_s expansion truncated at P_1; mean
cosine `f_1` is the parameter). General-anisotropic extension cited
to Mika 1961 but not pursued.

**Validity range**: c ≤ 1 + 1/(3 f_1) (Eq. 5). Outside this band the
discrete-mode count grows and the bi-orthogonality relations break.
For f_1 = 0 this reduces to "all c" (the isotropic limit). For
f_1 = 0.30 the upper bound is c ≤ 2.111.

# 2. Generalized dispersion + slab criticality

- Discrete EFs and continuum EFs: Eqs. 7–10.
- Dispersion (purely imaginary ν_0 for real c > 1): Eq. 11; quadratic
  in c: Eq. 12.
- Slab criticality condition (linearly anisotropic, **even modes** /
  symmetric flux): Eq. 46. The form is
  ±π/2 − arctan{ R sin[(d−z_0)/|ν_0|] + sin[(d+z_0)/|ν_0|] } /
                { R cos[(d−z_0)/|ν_0|] − cos[(d+z_0)/|ν_0|] }
  =          arctan{ … numerator with K_0,K_1,K_2 … } /
                { … denominator with K_0,K_1,K_2 … }
  where K_j = ∫₀¹ dν ν^j g_1(c,ν)/ν̄ · (cν/2) X²(−ν) Λ(∞) d(−νν̄) T(R,ν).
- T(R,μ) defined by Eq. 33; X(μ) by Eq. 40 (numerically friendly form);
  z_0 (extrapolated endpoint) by Eq. 42.

# 3. Sphere problem — the slab/sphere unification

Atalay reuses the slab machinery for the sphere via the standard
Case-Zweifel **antisymmetric trick**: the sphere flux corresponds to
the **odd-mode** slab solution. Specifically:
- ψ(x,μ) = −ψ(−x,−μ) (Eq. 47), giving a_{0+} = −a_{0−} and
  A(ν) = −A(−ν).
- T₁(R,μ) replaces T(R,μ) (sign flip in numerator, Eq. 50).
- L_j integrals replace K_j (Eq. 51).
- Sphere criticality: Eq. 54. Same machinery, single sign flip in two
  T-style functions and a sin↔cos shuffle. **Confirmed: this is the
  Siewert–Thomas sphere=slab+sign-flip pattern**.

This means **one implementation covers both geometries** with a
geometry flag.

# 4. Results extracted

- **Table 1**: ν̄(c) and z_0(c) for f_1 = 0.0 only, c ∈ [1.001, 8.10].
  (Surprising — anisotropic z_0(c, f_1) tables NOT given.)
- **Tables 2–5**: slab critical thicknesses 2d (mfp), even modes,
  for f_1 ∈ {0.0, 0.10, 0.20, 0.30}, R ∈ {0, 0.25, 0.50, 0.75, 0.99},
  c ∈ {1.01, 1.10, …, 2.00}. First three eigenvalues stacked per cell
  (the second/third often missing for c ≈ 1.01 or large f_1).
- **Tables 6–9**: c eigenvalues for fixed slab thickness 2d ∈
  {0.20, 2.0, 20.0}, f_1 ∈ {0.0, 0.10, 0.20, 0.30}, **even modes**.
- **Table 10**: c eigenvalues, **odd modes** (= sphere), f_1 = 0.10 only.
- **NO R = 1.0 column** — perfect reflection makes thickness drop out
  of the criticality condition (mentioned p.246).
- **NO sphere R-table at f_1 = 0.0, 0.20, 0.30** — only Table 10.

# 5. Mapping to Sood et al. 2003 cases

Sood naming convention `M-G-N(thickness)-A-B-G` where A=BC, B=anisotropy,
G=geometry. Atalay's data directly populates:

| Sood case prefix       | Atalay table         | Notes                              |
|------------------------|----------------------|------------------------------------|
| `Ua-1-0-SL`, `PUa-1-0-SL` (vacuum, isotropic) | Table 2 R=0      | Cross-check vs Sood [21,26]    |
| `Ua-1-1-SL` (vacuum, anisotropic)             | Tables 3-5 R=0   | NEW unlock for f_1>0           |
| `*-1-N(2d)-1-0-SL` (reflected, isotropic)     | Table 2 R>0      | NEW unlock                     |
| `*-1-N(2d)-1-1-SL` (reflected+anisotropic)    | Tables 3-5 R>0   | **NEW — cross-product unlock** |
| `Ua-1-0-SP`, `PUa-1-0-SP` (vacuum sphere, isotropic)  | Table 6 R=0  | Cross-check |
| `*-1-N-1-0-SP` (reflected sphere, isotropic)  | (only Table 10 covers sphere — f_1=0.10 only) | PARTIAL |
| `*-1-N-1-1-SP` (reflected+anisotropic sphere) | Table 10 only    | f_1=0.10 ONLY                  |

**Cross-product (reflected+anisotropic) cases ARE in Sood** (e.g.,
Engle-Mendelson "U-1-2-SL" with linear anisotropy + reflector).
Atalay populates the slab cross-product densely (Tables 3-5 ×
all R), but the **sphere cross-product is single-row** (only f_1=0.10).

# 6. Comparison with the other 3 papers

- **Burkart-Ishiguro-Siewert 1976**: F_N method, linearly anisotropic,
  vacuum-only slab+sphere. **Complementary**: F_N is
  collocation/numerical-quadrature-cleaner; Atalay is exact-Fredholm-
  first-order. Use B-I-S as cross-check on `*-1-1-SL` (R=0).
- **Dahl–Sjöstrand 1979**: Carlvik high-order spatial expansion,
  linearly anisotropic eigenvalue spectrum. Atalay says results
  "agree closely" (p.246) — useful as a third independent
  cross-check on Tables 6–10.
- **Neshat–Maiorino**: not in this PDF, no comment.

**Verdict — Atalay 1997 IS the primary source** for the
reflected+anisotropic combined matrix because:
1. It's the only paper that systematically tabulates f_1 ∈ {0, 0.1,
   0.2, 0.3} × R ∈ {0, 0.25, 0.5, 0.75, 0.99} × c grids for slab.
2. It's the only paper using a UNIFIED slab/sphere derivation.
3. The Fredholm framework is iterable to higher accuracy
   (p.246 — first order is published, higher orders feasible).

Use B-I-S, D-S, Kohut as **cross-check / triangulation** sources.

# 7. Implementation feasibility — fn_method/

**Critical finding**: Atalay is **NOT F_N method** — it is the
**Case singular-eigenfunction + Fredholm iteration** branch. It does
NOT extend `fn_method/core/` cleanly. The X-function (Eq. 26 / 40),
γ(μ) (Eq. 23), and N(ν) (Eq. 36) live in a different mathematical
machinery than the F_N matrix-collocation code path.

**Recommended architecture**:
```
case_method/                       # NEW package (NOT inside fn_method/)
├── core/
│   ├── x_function.py              # Eq. 40
│   ├── extrapolated_endpoint.py   # Eq. 42
│   └── dispersion.py              # Eqs. 11, 12
├── slab/                          # Eq. 46 (even), Eq. 39 (general)
├── sphere/                        # Eq. 54 (= odd-mode slab via Eq. 47)
└── anisotropy/
    └── linear.py                  # f_1 ≠ 0 corrections to all kernels
```

The slab/sphere split is by EVEN/ODD mode, NOT by physical geometry —
sphere reuses 100% of slab kernels with one BC sign flip + a redefined
T → T_1. This is the **same architectural insight as Siewert-Thomas
1986** for sphere via slab.

**SymPy walls**: X(μ) integral (Eq. 40) has a logarithmic singularity
at ν = μ — Atalay says "successive subdivision near right endpoint"
(p.246 footnote). Closed-form NOT available; needs mpmath quadrature.
The dispersion roots (Eq. 12) are clean quadratic in c — SymPy can
do those.

**Tractability**: implement in **two passes**:
1. **Pass 1**: f_1 = 0 (isotropic) cross-check vs Atalay 1995 + Sood
   isotropic table — exercises the X-function + extrapolated-endpoint
   + Fredholm machinery without anisotropy bookkeeping.
2. **Pass 2**: f_1 > 0 — adds the K_0/K_1/K_2 (and L_0/L_1/L_2)
   integrals on top of validated Pass-1 infrastructure.

# 8. Recommended Sood-case anchoring

| Case class                            | PRIMARY        | CROSS-CHECK              |
|---------------------------------------|----------------|--------------------------|
| Vacuum slab, isotropic, 1G            | Sood [21,26]   | Atalay R=0, Kaper 1974   |
| Vacuum slab, anisotropic              | **Atalay** R=0 | Burkart-Ishiguro-Siewert  |
| Reflected slab, isotropic             | **Atalay** f_1=0 | Garis-Sjöstrand 1994    |
| Reflected slab, anisotropic           | **Atalay**     | Dahl-Sjöstrand, Kohut    |
| Vacuum sphere, isotropic              | Sood + Siewert-Thomas | Atalay Table 6 R=0 |
| Vacuum sphere, anisotropic            | **Atalay** R=0 (Tables 7-10) | none clean    |
| Reflected sphere, isotropic           | **Atalay** Table 6 R>0 | none           |
| Reflected sphere, anisotropic         | **Atalay** Table 10 (f_1=0.1 only) | NO independent source — gap |

# 9. Errata watch

Spot-check on Eq. 32 vs Eq. 34: the substitution of T(R,μ) and
γ(μ) into Eq. 27 to produce Eq. 34 needs the κ-dependent factor
ν^2 d(ν_0 ν̄) + ν_0 ν̄ + ν ν̄ d(ν_0²) − ν_0² inside the brackets —
this matches between Eqs. 32 and 34 (verified by inspection).
Eq. 39 numerator/denominator structure is symmetric under the
K_j → K_j' reflection, consistent with the complex-conjugate
argument used to derive Eq. 46. No errata flagged in a quick scan.

The R = 1 limit makes the slab thickness drop out of Eq. 46
(p.246) — that's a feature of the physics (perfect reflector =
infinite medium), not a bug.

# 10. Suggested user action

Add Atalay 1997 to Zotero (DOI 10.1016/0149-1970(95)00094-1). Tag
with `case-method`, `reflected-bc`, `anisotropic-scattering`,
`slab+sphere-criticality`. The paper is short (24 pp.), notation
is internally consistent, and Tables 2–10 are essentially a
ready-made benchmark suite for the reflected+anisotropic Sood
gap.
