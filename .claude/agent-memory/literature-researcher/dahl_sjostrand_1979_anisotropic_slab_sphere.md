---
name: Dahl-Sjostrand 1979 — anisotropic slab+sphere eigenvalue spectrum
description: NSE 69, 114-125. Carlvik integral-equation Legendre-Galerkin spectral solver for c-eigenvalue spectrum of bare slab/sphere with linearly anisotropic scattering. Complete Tables I/II up to 11 eigenvalues × {μ̄=0, 0.10, 0.15, 0.20, 0.25, 0.30} × d∈{0.2..20}. Verdict: PRIMARY for both anisotropic slab and sphere; cross-check baseline against BIS-76 (slab only, F_N).
type: reference
---

# Dahl & Sjöstrand (1979), *NSE* **69**, 114-125 — DOI: 10.13182/NSE69-114

## 1. Method classification — NOT F_N

It is a **Legendre-Galerkin spectral expansion of Carlvik's integral
equation**. Direct Fredholm-like discretization on Legendre polynomials
of the spatial variable, NOT singular-eigenfunction (Case) and NOT
F_N (boundary-collocation). Structurally INDEPENDENT of every F_N
result in Sood — kills any cross-method correlation risk.

**Eq. (1) Carlvik integral equation** (verbatim, p. 118):

φ(x) = (c·a/2) { ∫_{-1}^{+1} φ(y) E_1(a|x-y|) dy
              − 3μ̄(c-1) [ ∫_{-1}^{+1} φ(y) E_3(a|x-y|) dy
                         − (E_3(a|1-x|) + (-1)^q E_3(a|1+x|))/2
                           · ∫_{-1}^{+1} φ(y) y^q dy ] }

with a = half-thickness or radius in mfp, q=0 (slab) / q=1 (sphere),
φ = flux density (slab) or r·φ (sphere). Vacuum BC.

**Eq. (2) Legendre expansion**: φ(x) = Σ F_n (2n+1) P_n(x) — even n
for slab (fundamental even), odd n for sphere (r·φ is odd → P_1, P_3,
… P_17).

**Eq. (3) Matrix eigenvalue problem**:

[A − 3μ̄(c−1) B] F = (1/(c·d)) F .

A_{m,n}, B_{m,n} closed-form recurrences (Carlvik 1968, NSE 31, 295;
**erratum noted: sign of last term in Carlvik Eq. 4b is misprinted**).

**Eq. (4) Linearization to standard eigenproblem** (the trick that
makes it cheap):

| G  H | (F)        (F)
|      |     = (1/c)
| I  0 | (K)        (K)

with G = d(A + 3μ̄B), H = −3μ̄·d·B, K = c·F. Doubles matrix size but
removes the c-dependence inside the brackets → ALL eigenvalues found
in one standard call. Carlvik's iterative scheme is replaced.

## 2. Mapping to Sood LA-13511

**Anisotropic order**: P_1 only (μ̄·c-1 in front of E_3 term). Does
NOT cover P_2.

**Materials**: c-eigenvalue formulation is material-AGNOSTIC. To map
to Sood `*-1-1-SL/SP` cases, look up that case's (c, μ̄, σ_t·R)
triple, find the d=2σ_t·R row in the table, and read c off Eigenvalue
1. Geometry/material decoupling means **every** Sood `*-1-1-SL/SP`
benchmark with μ̄ ∈ {0.0, 0.10, 0.15, 0.20, 0.25, 0.30} is
covered — that is the full P_1 family.

**P_2 coverage**: NONE. Falls back to BIS-76 or Atalay-1997 for the
`*-1-2-SL/SP` cases. The paper's framework is structurally limited
to linearly anisotropic; Eq. (1) only carries one E_3 anisotropic
correction term.

## 3. Method specifics — what it actually delivers

- **Spectrum, not just dominant**: Tables I (sphere) and II (slab)
  give **up to 11 eigenvalues** per (d, μ̄). Real eigenvalues for
  μ̄=0 (isotropic). For μ̄ ≥ 0.15, complex-conjugate pairs appear at
  higher modes — flagged with ± and i (e.g. d=0.2, μ̄=0.30:
  eigenvalue 3 = 24.9156i, eigenvalue 4 = 66.5149 ± something).
  **Higher harmonic eigenvalue spectrum is the paper's distinctive
  contribution** — BIS-76 gives only the dominant mode.
- **Numerical pathology**: complex eigenvalues at high anisotropy.
  Real part of complex eigenvalues DECREASES as μ̄ increases — the
  pairs migrate down in c. With increasing μ̄, complex eigenvalues
  extend to LOWER mode index AND LARGER bodies (Figs. 3-6 show this
  clearly). For ORPHEUS: a fundamental-mode k_eff verifier doesn't
  hit this, but a multi-mode α-eigenvalue verifier does.
- **Convergence**: 9 Legendre terms (P_0..P_16 even slab, P_1..P_17
  odd sphere). 10-term run agrees to "the number of digits given
  in Tables I and II" → max 1-unit error in last digit, except
  underlined entries (10-unit error). Excellent precision.
- **Self-consistency**: agreement with Hartman 1975 (isotropic) is
  excellent; with Kerner et al. 1967 also good. Anisotropic
  fundamental matches Kschwendt 1971 and Hembd-Kschwendt 1972.

## 4. Implementation feasibility

- **Pillar 2 / structurally independent from Sood F_N truth set**:
  YES — Galerkin on integral form, not F_N. Safe redundancy.
- **SymPy walls**: A_{m,n}, B_{m,n} closed-form recurrences from
  Carlvik 1968 — symbolic-OK at low rank. E_1, E_3 = exponential
  integrals: SymPy has `expint(n,x)`, but for d ≥ 5 needs mpmath
  fallback (numerical conditioning of high-order recurrences, per
  Carlvik 1968 erratum).
- **Generalized eigenproblem**: Eq. (4) is a standard non-symmetric
  GEP — `scipy.linalg.eig` covers it; complex eigenvalues land
  natively.

## 5. Architectural verdict — PRIMARY for sphere AND slab P_1, paired with BIS-76 as cross-check

- **PRIMARY (P_1, both geometries)**: Dahl-Sjöstrand. Why:
  (a) covers BOTH slab and sphere in one paper / one method;
  (b) full eigenvalue spectrum (≥6 eigenvalues per case, often 11);
  (c) tabulated to 7 significant digits;
  (d) cleanest math (one transcendental → one linear GEP).
- **CROSS-CHECK (P_1 slab)**: BIS-76. Structurally different (F_N
  collocation vs Galerkin) — independent verification.
- **CROSS-CHECK (P_1 sphere)**: NONE in this corpus. Sphere
  validation rests on Dahl-Sjöstrand alone unless a separate sphere
  F_N reference (e.g. Sanchez 1986 NSE for transport in a homogeneous
  sphere with linearly anisotropic scattering) is added.
- **P_2**: Dahl-Sjöstrand does NOT cover this — cleanly falls back
  to BIS-76 / Atalay-1997.

## 6. Errata noted in the paper itself

Authors flag: **Carlvik 1968 NSE 31, 295, Eq. (4b) — sign of last
term is misprinted**. Use Dahl-Sjöstrand recurrences as the corrected
master, not Carlvik 1968 directly.
