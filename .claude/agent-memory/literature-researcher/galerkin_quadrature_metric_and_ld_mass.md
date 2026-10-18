---
name: galerkin-quadrature-metric-and-ld-mass
description: Galerkin quadrature GQ1/GQ2/GQ3 (Morel 1989, Morel+ 2017, Shands+ 2025) and which Gram the moment axis carries (Sanchez-Ragusa 2011 W_gal = DᵗD); curvilinear LD modal mass matrix V,W,X; nodal row-sum = Palmer-Adams FL; weighted cellular sheaf (Hansen-Ghrist 2019)
metadata:
  type: project
---

Extracted 2026-10-18 for the numerics-untangling frame-bridge question. Memo:
`scratch/numerics_untangling/frame_bridge/literature.md`.

**Now LOCAL** (copied from the NSE zip archive; OCR got 429, so read via `pdftotext`; re-run
the OCR later):
- Sanchez-Ragusa(2011), NSE 169:133, `10.13182/NSE10-31`;
- Morel-et-al(2017), NSE 185:325, `10.1080/00295639.2016.1272383`;
- Shands-et-al(2025), NSE 199(5):854, `10.1080/00295639.2024.2385220`;
- Hansen-Ghrist(2019), arXiv 1808.01513;
- Arnold-Brezzi-Cockburn-Marini(2002), SINUM 39:1749 (author copy from umn.edu).

**Galerkin quadrature**
- **GQ1**: Morel 1989 Eqs. 24c and 27. D is the continuum moment of the interpolant;
  M = D⁻¹.
- **GQ2**: Morel 2017 Eq. 32. D is the standard P(μ)w; M = D⁻¹. It gives exact P_{N−1}
  moment equations with a closure.
- **GQ3**: Shands 2025 Eqs. 35–40. Legendre vectors are Gram-Schmidt-orthogonalised under
  the quadrature, so DM = I with no inversion.
- **Sanchez-Ragusa Eq. 21**: the moments carry the continuum (orthonormal spherical
  harmonic) Gram. The ordinate scalar product is W_gal = DᵗD (p. 138).
- `[M]` W_gal equals the exact consistent angular mass ∫B_mB_n. It is diagonal only for
  Gauss.
- GQ1 scattering is NOT diag(w)-self-adjoint for double-Gauss.
- No published adjoint-consistency analysis of GQ was found. That gap is real as far as I
  searched.

**LD mass**: Adams-Martin (A.4a–c) are verified on the rendered page. The sphere
off-diagonal entry is r_kΔ²/3, so the matrix is never diagonal. (A.11a): the unknown f_0k
is the Δr-mean, and the volume mean is f_0k + (W/V)f^x.

**Palmer-Adams 1993 FL** T_k and R_k = the ROW SUMS of the nodal consistent matrices
(SymPy-verified). This is why the curvilinear literature is nodal: lumping lives in the
Lagrange basis.

**Hansen-Ghrist 2019**: §3.1 weighted sheaf (an inner product per stalk); §3.5 discrete
vector bundle (all restriction maps invertible). Related: [[issue-158-linear-discontinuous-cell-update]].
