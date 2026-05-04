# Issue #113 close-out — Peierls slab K-matrix failed-scheme forensic

**Status:** CLOSED. The naive GL-collocation K-matrix assembly (one-point cross-panel rule + singularity-subtraction on the diagonal) has been retired. Production path is now the unified basis-aware adaptive quadrature in `orpheus.derivations.peierls_geometry.solve_peierls_1g` / `solve_peierls_mg` (the `_pg.solve_peierls_*` family).

**Hard evidence pins:** ERR-027 and ERR-028 in `tests/l0_error_catalog.md`.

**Date:** 2026-04-30
**Cleanup commit:** <!-- COMMIT-HASH -->
**Source removed from docs:** `docs/theory/collision_probability.rst`, lines ~2324–2372 (the user-quoted range 2324–2347 covers the diagnostic; the immediately-following fix narrative through line ~2372 is included here for self-containment).

---

## Surviving Sphinx anchors

The forensic narrative below is removed from `collision_probability.rst`. The surviving theory-page presentation of the production path lives in `docs/theory/peierls_nystrom.rst`:

- `:eq:\`peierls-slab-polar\`` (line ~536) — the slab polar-form equation that the unified K-matrix discretizes.
- `:eq:\`peierls-unified\`` (line ~2821) — the geometry-agnostic Peierls operator.
- `:eq:\`peierls-vacuum-bc-slab\`` (line ~3426) — vacuum-BC verification milestone (machine-precision K-matrix gate).
- `:eq:\`peierls-white-bc-slab\`` (line ~3550) — White-BC slab analytical flux (Wigner–Seitz identity).
- Subsection "The slab polar-form equation (recap)" (`peierls_nystrom.rst` line ~528) — the surviving K-matrix construction reference.
- Subsection "Why the outer μ-integral is stiff at μ = 0" (`peierls_nystrom.rst` line ~571) — surviving rationale for adaptive quadrature.

---

## Failed-scheme forensic (verbatim, relocated from `collision_probability.rst`)

### Setup — composite Gauss–Legendre panels

The slab $[0, L]$ is divided into $n_{\rm panels}$ panels per material region. Within each panel, $p$ Gauss–Legendre nodes and weights are computed. The total number of quadrature nodes is $N = n_{\rm regions} \times n_{\rm panels} \times p$.

### Singularity treatment — classical decomposition

The $E_1$ kernel has a logarithmic singularity at $x = x'$:

```math
E_1(z) = \bigl[-\ln z - \gamma\bigr] + R(z),
\qquad R(z) \equiv E_1(z) + \ln z + \gamma,\quad R(0) = 0.
```

The remainder $R(z)$ is a smooth (analytic) function that vanishes at the origin. This decomposition motivates the classical **singularity-subtraction** approach (used in the original implementation), in which diagonal panels split the kernel into a smooth $R$ part (handled by GL weights) and a $-\ln|x_i - x'|$ part (handled by product-integration weights). Off-diagonal panels would then use standard GL weights on the smooth $E_1$ integrand.

### Why this classical scheme fails in practice (issue #113, ERR-027/028)

The scheme rests on two assumptions that turn out to be insufficient:

1. *Off-diagonal GL is exact for smooth integrands.* This is true for polynomial integrands of degree $\le 2p-1$, but $E_1$ is transcendental with a near-log spike at $x'\to x_i$. One-point collocation (`E_1(τ_ij)·w_j`) gives ~1% error even for panel pairs with moderate optical separation, worst at panel-boundary neighbours where the log spike sits just outside the source panel.

2. *R is smooth across the diagonal panel.* True in $\tau$ but $R(\Sigma_t |x_i - x'|)$ has a derivative kink in $x'$ at $x'=x_i$ (from the absolute value). GL cannot integrate across interior derivative discontinuities. This adds ~1% error on every diagonal panel.

Both bugs were invisible to row-sum conservation tests because $\sum_j L_j(x')=1$ kills the kink in the summed integrand even when each basis-individual $K[i,j]$ is wrong — hence the "passes at mode 0, fails at mode n>0" signature observed in the rank-N investigation.

### Diagnostic signature

- Cross-panel $K[i,j]$ entries: 0.5–1.5% relative error at $n_{\rm panels}=2$, $p=4$; refining panels reduces this only at $\mathcal{O}(h)$ rate (non-spectral). Worst case: $x_i$ within a small optical distance of the source-panel boundary.
- Diagonal-panel entries (~40% of all $K$ entries for $p=4$, $n_{\rm panels}=8$): ~1% relative error; element-wise $K[i,i]$ deviates by ~1.2e-2 at $n=2$, $p=4$.
- Downstream: slab k-eff tests showed ~0.4% tie-point offsets in the Sanchez R≈1.98 case (Sanchez1982).
- Row-sum $K \cdot [1, 1, \dots]^\top$ stays exact to $\mathcal{O}(10^{-15})$ — partition-of-unity of the Lagrange basis cancels the kink in the summed integrand. Every existing row-sum test was blind. The "passes at mode 0, fails at mode n>0" signature appeared only when basis-individual modes were probed (rank-N investigation).

### The fix — unified basis-aware adaptive assembly

Production implementation: `orpheus.derivations.peierls_slab._basis_kernel_weights`, called from the unified `_pg.solve_peierls_*` path. Every $K[i, j]$ is computed directly as

```math
K[i, j] \;=\; \tfrac{1}{2} \int_{p_a}^{p_b}
              E_1\!\bigl(\tau(x_i, x')\bigr)\,L_j(x')\,\mathrm{d}x'
```

via adaptive `mpmath.quad`, with the subdivision hint $[p_a, x_i, p_b]$ when $x_i$ lies inside the source panel (same-panel case). This single code path:

- handles the integrable log singularity at $x'=x_i$ via the subdivision hint (mpmath resolves it to machine precision);
- handles the derivative kink of $R(\tau)$ in $x'$ — GL on each smooth half of the subdivided panel converges spectrally;
- eliminates the off-diagonal-panel quadrature error — adaptive refinement resolves $E_1$'s non-polynomial structure on arbitrary panel pairs without relying on near-log assumptions.

The implementation exactly mirrors the adaptive reference `orpheus.derivations.peierls_reference.slab_K_vol_element`, so the production code and the reference agree to machine $\mathrm{dps}$ by construction.

### Why adaptive quadrature over the alternatives

Four strategies exist for the Peierls log singularity (Atkinson1997):

1. *Graded meshes* (cluster GL nodes near the diagonal): algebraic convergence only; many nodes needed for high accuracy.
2. *IMT transformation* (Iri–Moriguti–Takahashi double-exponential): truncated here — see the production rationale on `peierls_nystrom.rst`.

(Full enumeration of strategies 3–4 was preserved in the surviving "Why the outer μ-integral is stiff at μ = 0" subsection of `peierls_nystrom.rst`.)

---

## Hard evidence pins (verbatim from `tests/l0_error_catalog.md`)

### ERR-027 — Peierls slab K-matrix: naive GL collocation for cross-panel entries

- **Failure mode:** #3 Missing factor — missing quadrature resolution (one-point rule where adaptive is required)
- **Date:** 2026-04-19
- **Solver:** CP / Peierls (`orpheus.derivations.peierls_slab._build_kernel_matrix`).
- **L1 test that catches it:** `tests/derivations/test_peierls_reference.py::TestSlabKMatrixElementwiseVsReference::test_cross_panel_boundary_neighbour_elementwise` — element-wise `K[4, 3]` at n_panels=2, p=4 vs the adaptive `slab_K_vol_element` reference at 1e-10.
- **Lesson:** "Row-sum conservation" tests are systematically blind to basis-individual quadrature errors that happen to sum to zero under partition-of-unity of the Lagrange basis. Element-wise $K[i, j]$ verification against an adaptive reference is the only reliable L0/L1 gate for Nyström kernel assembly.

### ERR-028 — Peierls slab K-matrix: GL collocation of remainder R(τ) has unresolved kink at x'=x_i

- **Failure mode:** #3 Missing factor — missing subdivision hint
- **Date:** 2026-04-19
- **Solver:** CP / Peierls (`orpheus.derivations.peierls_slab._build_kernel_matrix`).
- **L1 test that catches it:** `tests/derivations/test_peierls_reference.py::TestSlabKMatrixElementwiseVsReference::test_small_case_elementwise_agreement` — element-wise $K[i, j]$ including diagonal entries vs the adaptive `slab_K_vol_element` reference at 1e-10.
- **Lesson:** Singularity subtraction splits a kernel into "smooth" and "singular" parts — but "smooth in τ" is not the same as "smooth in x'". A change of variables that folds the singularity into a kink instead of removing it merely trades one unresolved feature for another. The adaptive-with-hint approach is more robust and unifies all cases.

Related downstream catalog entry: ERR-029 (Peierls curvilinear K-matrix: ρ/ω integration does not subdivide at ray-panel crossings or tangent angles) — same family of "missing subdivision hint" failure mode in the curvilinear path.

---

## Production decision

- The naive scheme (one-point cross-panel collocation + singularity-subtraction on the diagonal) is **retired**.
- The unified adaptive-quadrature K-matrix assembly (`_basis_kernel_weights` + `mpmath.quad` with subdivision hint) is the **single production path** for slab, cylinder, and sphere; geometry-specific differences are factored behind `CurvilinearGeometry` primitives.
- Element-wise verification (not row-sum) is the canonical L0/L1 gate going forward — pinned by the two `TestSlabKMatrixElementwiseVsReference` tests above.

Citations: Sanchez1982, Stamm1983 (cited in the surviving `peierls_nystrom.rst` for the unified scheme); Atkinson1997 (quadrature-strategy comparison).
