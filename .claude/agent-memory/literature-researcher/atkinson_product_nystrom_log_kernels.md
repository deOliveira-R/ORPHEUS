---
name: Atkinson product-Nyström for log-singular kernels
description: Atkinson 1997 book §4.2 + Atkinson 1972 Numer. Math. paper, both LOCAL — full method spec for product-Simpson/trapezoidal on log|t-s| Fredholm kernels, with graded-mesh and de Hoog-Weiss superconvergence
type: reference
---

## Where it lives

Both Atkinson sources are local in `scratch/literature/`:

- `Atkinson - The numerical Solution of Integral Equations of the Second Kind/`
  — 16 PDFs (one per chapter). Chapter 4 = `07.0_pp_100_156_The_Nystroem_method.pdf`
  — §4.2 "Product integration methods" pp. 116–135 is the canonical method spec.
- `Atkinson(1972)The numerical solution of Fredholm integral equations of the
  second kind with singular kernels.pdf` — 12 pp., the practical
  hybrid product/standard Simpson algorithm with adaptive near-/far-field
  splitting at threshold $\delta$.

## What it gives

For $\lambda x(t) - \int K(t,s) x(s) ds = y(t)$ with $K = L\,H$
(L smooth, H = log|t-s| or |t-s|^{γ-1}):

1. **Method**: piecewise-polynomial interpolation of $L\cdot x$ in $s$;
   integrate the singular factor analytically against the basis
   (Eq. 4.2.63–4.2.67 product-trapezoidal; §4.2.1 product-Simpson).
2. **Weight tables**: ψ_0(k), ψ_1(k) closed-form integrals (Eq. 4.2.74)
   give all $w_j(t_i)$ in $O(N)$ setup cost. Implementation gotcha at
   p. 119: large-|k| cancellation requires asymptotic series.
3. **Convergence**: Theorem 4.2.1 + de Hoog-Weiss (Eq. 4.2.83) gives
   $O(h^4 \log h)$ for product Simpson on log|t-s| kernels — IF
   $x \in C^4$.
4. **Reality check (Schneider Theorem 4.2.2, p. 126)**: solutions of
   $(\lambda - \mathcal K)x = y$ with log kernel have endpoint
   singularities $(t-a)^\gamma$, killing high-order rates on uniform
   meshes.
5. **Graded-mesh repair (Theorem 4.2.3, p. 132)**: $t_j = (j/n)^q$
   with $q > m+1$ for log kernels, $q \ge (m+1)/\gamma$ for
   $|t-s|^{\gamma-1}$. Restores $O(n^{-(m+1)})$.

## Bickley-Naylor connection (NOT in Atkinson — derived in memo)

$\mathrm{Ki}_1(\tau) = \pi/2 + \tau(\log\tau + \gamma_E - 1 - \log 2) + O(\tau^2)$,
so the slab Peierls kernel $\frac12 \mathrm{Ki}_1(\Sigma|x-y|)$ falls
under Atkinson's (4.2.80) decomposition with $L_2 = (\Sigma/2)|x-y|$,
$H_2 = \log(\Sigma|x-y|)$ + smooth Taylor remainder treated as
$L_1 H_1$ with $H_1 \equiv 1$. The sphere Peierls kernel
$|r-r'|^{-1}\mathrm{Ki}_1(\Sigma|r-r'|)$ also reduces to $\log|r-r'|$
behavior because $\mathrm{Ki}_1(\tau) \sim \tau\log\tau$ near 0.

## Reference network (Atkinson book pp. 154-156)

- Young 1954 [ref 585] — original product-integration idea.
- de Hoog & Weiss 1973 [ref 164] — superconvergence for log/algebraic.
- Schneider 1979 [ref 493, 494] — solution regularity + graded mesh.
- Rice 1969 [ref 461] — $(j/n)^q$ graded mesh.
- Graham 1981 [ref 227] — Hölder-space extensions.

## Memo location

`scratch/derivations/peierls_log_singular/atkinson_product_nystrom.md`
— full method spec with notation bridge, equation-by-equation, plus
implementation pseudocode for `flux_reconstruction.py`.

## Use case

When numerics-investigator picks up Wave 2-A "Path A.i flux" stall
(F_N moments converged to 1e-5 but reconstructed scalar flux stalls
at 1-3%): Atkinson-graded product Simpson with q=4 is the textbook fix.
