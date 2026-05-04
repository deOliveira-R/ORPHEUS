# Issue #133 close-out — Phase 5 continuous-µ retreat (closed wontfix)

**Date:** 2026-04-30
**Doc cleanup commit:** <!-- COMMIT-HASH -->
**Source relocated from:**
- `docs/theory/peierls_nystrom.rst` lines 7017–7196 (Phase 5a + Round 1 (3 fronts) + Round 2 (PRIMARY+BACKUP) + Round 3 (PRIMARY+SECONDARY)).
- `docs/theory/peierls_nystrom.rst` lines 9665–9687 (§22.7 visibility-cone "Provenance" historical detail).

The §22.9 quadrature-rollout content for #133 specifically (what landed in the chord_quadrature recipe) is the responsibility of a separate close-out comment for Issues #133/#134/#135/#136.

---

## Summary

- The continuous-µ reformulation of $K_{\rm bc}$ for the sphere specular boundary appeared to be the structural fix for the matrix-Galerkin divergence at high rank $N$. **It is not.** The kernel is **hypersingular** at the discrete Nyström diagonal and cannot be discretised by any standard quadrature technique.
- 6 rounds of investigation (Phase 5a, Round 1 with 3 fronts, Round 2 with PRIMARY+BACKUP, Round 3 with PRIMARY+SECONDARY) independently confirmed the obstruction. The matrix-Galerkin form $(I - T R)^{-1}$ and the continuous-µ integral form are **different integral operators**, not different discretisations of the same operator. The matrix-Galerkin's mode-mixing via the basis projection is the load-bearing structure giving Phase 4 its operational behaviour — Nyström sampling cannot reproduce it.
- **Production decision (preserved):** `compute_K_bc_specular_continuous_mu_sphere` is shipped as a reference implementation with SymPy 4/4 verification, but the closure `boundary="specular_continuous_mu"` is registered in dispatch and **raises `NotImplementedError`** with a guidance message pointing to `closure="specular_multibounce"` (Phase 4 matrix-Galerkin form) for production use. The earlier Phase 4 docstring reference to "Phase 5 — proper fix" is hereby withdrawn: there is no proper fix in the continuous-µ discretisation framework.
- **Useful side-finding (promote-worthy):** the chord/visibility-cone substitution $u^2 = (\mu^2 - \mu_{\min}^2)/(1 - \mu_{\min}^2)$ gives **machine-precision off-diagonal Q-convergence** (1e-9 at Q=128) for any µ-resolved per-pair integral with a single-endpoint visibility-cone singularity. Portable to any geometry. This substitution survives in the production `chord_quadrature` recipe (see §22.9 in the theory page).

Surviving doc anchors:
- `:label:` `peierls-phase5-sanchez-A6` — Sanchez 1986 Eq. (A6), the textbook continuous-µ form, retained for documentation.
- `:ref:` `peierls-phase5-retreat` — section anchor (the section title that survives as a stub pointing here).
- `[SanchezTTSP1986]` citation entry preserved in the references.

<details><summary>Round-by-round full investigation record (verbatim from peierls_nystrom.rst:7017–7196)</summary>

### Phase 5 — continuous-µ reformulation (RESEARCH RECORD, ABANDONED for production wiring after 3 rounds of investigation, 2026-04-28)

The continuous-µ reformulation appeared to be the structural fix for the matrix-Galerkin divergence at high rank $N$. It is not — the kernel is **hypersingular** at the discrete Nyström diagonal and cannot be discretised by any standard quadrature technique. The matrix-Galerkin $(I - T R)^{-1}$ mode-mixing absorbs the singularity via basis projection and is **essential**, not a numerical artefact.

### Round-by-round summary

#### Phase 5a (Sanchez 1986 reference implementation)

The textbook Eq. (A6) form is

```math
g_h(\rho' \to \rho) = 2 \int_{\mu_0}^{1} T(\mu_-)
    \mu_*^{-1} \cosh(\rho\mu) \cosh(\rho'\mu_*) e^{-2a\mu_-} d\mu
```

SymPy 4/4 verifications PASS. Reference implementation `compute_K_bc_specular_continuous_mu_sphere` in `orpheus/derivations/peierls_geometry.py`. Smoke test failed — kernel magnitude vs `closure="white_hebert"` ratios spread 4 orders of magnitude. Two open blockers: Sanchez↔ORPHEUS Jacobian conversion + diagonal singularity.

#### Round 1 (3 parallel fronts)

- **Front A (Jacobian conversion via rank-1 cross-check) FALSIFIED** — no separable conversion exists; Phase 5 kernel is rank-(>1), Hebert is rank-1, and Phase 5 itself doesn't converge in `n_quad`.
- **Front B (singularity subtraction) PARTIAL** — closed-form leading orders identified ($s = 2/[\mu(e^{\tau_0}-1)]$ interior, $s = 1/(a\mu^2)$ surface) but the analytical add-back is divergent.
- **Front C (ORPHEUS-native reformulation bypassing Sanchez cosh forms) FAILED** — built the $F_{\rm out} / G_{\rm in}$ µ-resolved primitives correctly (consistency passed vs Phase 4) but $k_{\rm eff}$ oscillates wildly with $Q$ because $T(\mu) \sim 1/(\sigma\cdot 2R\mu)$ at $\mu \to 0$ is the same singularity Phase 4 has at high $N$.

#### Round 2 (PRIMARY M2 + BACKUP symbolic limit)

- **PRIMARY M2 (bounce-resolved expansion)** — confirmed $w(\mu) = 1$ is the correct M1 weight, geometric series converges with ratio 0.25/bounce ($K_{\max}=10$ saturates), but the diagonal singularity persists at $K_{\max}=0$ (NO multi-bounce). Smoking gun: the singularity is in the $F_{\rm out} \cdot G_{\rm in}$ primitive structure, NOT in $T(\mu)$. Source: the $1/(\cos(\omega_i)\cos(\omega_j))$ Jacobian at $r_i = r_j = r$ has both cosines vanishing at the SAME $\mu_{\min}(r) = \sqrt{1 - (r/R)^2}$, yielding non-integrable $1/(\mu^2 - \mu_{\min}^2)$ on the visibility cone.
- **BACKUP (independent symbolic N→∞ limit)** derived the closed-form multi-bounce factor:

  ```math
  f_\infty(\mu) = \frac{1}{2} \frac{\mu}{1 - e^{-\sigma\,2R\mu}}
  ```

  bounded at $\mu = 0$ (limit $= 1/(4a)$); $(1/2)$ from $R = (1/2) M^{-1}$, $\mu$ from the µ-weighted basis Gram measure. Bit-exact via Bose-Einstein polylog: $K_\infty^{\rm half} = (1/(8a^2)) [\pi^2/6 - \mathrm{Li}_2(e^{-2a}) + 2a \ln(1 - e^{-2a})]$. The cross-domain attacker's M1 sketch was wrong by a factor of 2 (missing the $(1/2)$); Sanchez's Eq. (A6) is in a DIFFERENT normalisation (integral-equation Green's function vs Nyström $K_{ij}$).

#### Round 3 (PRIMARY adaptive quadrature + SECONDARY Galerkin)

6 approaches tried (per-pair half-M1, chord substitution $s^2 = \mu^2 - \mu_{\min}^2$, adaptive Gauss-Kronrod, Galerkin diagonal cell-average, full Galerkin double-integration $\int\int L_i(r) L_j(r') K(r, r') dr dr'$, alternative conventions). **All FAIL.** Smoke test $k_{\rm eff}$ errors −34 % to −50 %, monotone divergence in $n_{\rm quad}$. The full Galerkin double-integration over Lagrange basis (the standard BEM/IGA fix for integrable singular kernels) gives the same hypersingular log-divergence as the Nyström sampling — confirming the singularity is genuinely non-Nyström-discretisable.

### Structural conclusion

Independently confirmed by 3 numerics-investigators (Round 2 PRIMARY, Round 3 PRIMARY, Round 3 SECONDARY): the continuous-µ $K_{\rm bc}$ kernel is **hypersingular** (Hadamard finite-part / Cauchy principal-value type). The matrix-Galerkin form $(I - T R)^{-1}$ and the continuous-µ integral form are **different integral operators**, not different discretisations of the same operator. The matrix-Galerkin's mode-mixing via the basis projection is the load-bearing structure that gives Phase 4 its operational behaviour — Nyström sampling cannot reproduce it.

### Production verdict

`closure="specular_multibounce"` (Phase 4 matrix-Galerkin form) is the **permanent production path** for multi-bounce specular at all three geometries. Within its envelope ($N \in \{1, 2, 3\}$ for thin sphere/cyl, any $N$ for slab) it gives Hébert-quality $k_{\rm eff}$. The earlier Phase 4 docstring reference to "Phase 5 — proper fix" is hereby withdrawn: there is no proper fix in the continuous-µ discretisation framework.

### Useful side-finding (promote-worthy)

The chord/visibility-cone substitution $u^2 = (\mu^2 - \mu_{\min}^2)/(1 - \mu_{\min}^2)$ gives **machine-precision off-diagonal Q-convergence** (1e-9 at $Q=128$) for any µ-resolved per-pair integral with a single-endpoint visibility-cone singularity. Portable to any geometry.

### Sanchez 1986 reference theory (retained for documentation)

The Phase 5a reference implementation `compute_K_bc_specular_continuous_mu_sphere` produces Sanchez Eq. (A6) verbatim and the SymPy 4/4 verifications in `derivations/peierls_specular_continuous_mu.py` confirm the algebraic equivalence with the cross-domain attacker's M1 sketch (after the factor-of-µ correction). The closure is registered in the dispatch but raises `NotImplementedError` with a guidance message pointing to `closure="specular_multibounce"` for production use.

Below is the original Sanchez 1986 textbook derivation (kept as reference even though it does not yield a Nyström-compatible discretisation):

```math
g_h(\rho' \to \rho) \;=\;
     2\!\int_{\mu_0}^{1}\!T(\mu_-)\,\mu_*^{-1}\,
     \cosh(\rho\mu)\,\cosh(\rho'\mu_*)\,e^{-2a\mu_-}\,\mathrm d\mu
```
(label: `peierls-phase5-sanchez-A6`)

with $a = \Sigma_t R$, $T(\mu_-) = 1/(1 - e^{-2a\mu_-})$ the continuous multi-bounce factor, $\mu_-(\mu) = \sqrt{a^2 - \rho^2(1-\mu^2)}/a$ and $\mu_*(\mu) = \sqrt{\rho'^2 - \rho^2(1-\mu^2)}/\rho'$ dimensionless cosines, and chord-visibility cutoff $\mu_0 = \sqrt{\max(0, 1 - (\rho'/\rho)^2)}$. This bypasses the matrix-Galerkin projection entirely — the multi-bounce factor $T(\mu_-)$ is a multiplication operator, not the inverse of one, so there is no operator-norm divergence.

### Status of Phase 5a in ORPHEUS

- **Reference implementation shipped** as `compute_K_bc_specular_continuous_mu_sphere`. Calls Sanchez Eq. (A6) verbatim with $\alpha = 1, \beta = 0, \omega_1 = 0$ (perfect specular, isotropic scattering, homogeneous sphere).
- **SymPy verification shipped** at `derivations/peierls_specular_continuous_mu.py` — 4/4 checks PASS, including the µ-weight equivalence with the cross-domain attacker's M1 sketch (the M1 sketch had a factor-of-µ bug that SymPy V2 caught).
- **Sanchez 1986 literature memo** at `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md`.
- **Production wiring DEFERRED to Phase 5+**. The closure `boundary="specular_continuous_mu"` is registered in the dispatch but raises `NotImplementedError`. Two open blockers:

  1. **Diagonal singularity**: SymPy V4 surfaced that the integrand at $\rho' = \rho$ has a $1/\mu^2$ non-integrable singularity at the surface diagonal $\rho = a$ and a $1/\mu$ (logarithmic) singularity at interior diagonals. Sanchez does not specify a numerical method. ORPHEUS Nyström sampling at quadrature points hits the singularity directly; production wiring needs adaptive µ-quadrature, singularity subtraction, or a change-of-variables.
  2. **Sanchez ↔ ORPHEUS $K_{ij}$ Jacobian conversion**: Sanchez's $g_h(\rho' \to \rho)$ has the surface-area Jacobian $4\pi\rho'^2$ baked in via Eq. (2). ORPHEUS's discrete Nyström convention uses explicit `rv = 4π r²` radial-volume weights. The conversion factor is not closed-form in Sanchez 1986; deriving it (likely via a rank-1 cross-check against `boundary="white_hebert"`) is Phase 5+ work.

For multi-bounce specular today, use `boundary="specular_multibounce"` (Phase 4 matrix-Galerkin form). At rank-1 it reduces algebraically to `boundary="white_hebert"` for sphere/cyl, which has the same Hébert-quality accuracy as the Phase 5 closure would in the continuous limit.

### Reference

[SanchezTTSP1986] R. Sanchez, "Integral form of the equation of transfer for a homogeneous sphere with linearly anisotropic scattering," *Transport Theory & Statistical Physics*, vol. 14, pp. 333–343 (1986). DOI: 10.1080/00411458608210456.

</details>

---

## §22.7 visibility-cone substitution — Provenance (historical detail)

<details><summary>§22.7 Provenance subsection (verbatim from peierls_nystrom.rst:9665–9687)</summary>

The substitution surfaced during the Phase 5 Round 3 SECONDARY mission of the Peierls specular-BC investigation (commit `4dc03cf`, retreat note in the §"Phase 5 retreat" subsection above). The original diagnostic at `derivations/diagnostics/diag_phase5_round3_visibility_cone_quad.py` used a *linear-shift* variant $\mu = \mu_{\rm vis} + (1 - \mu_{\rm vis})\,u^{2}$ which works locally near the singular endpoint but does **not** give the global $\sqrt{y^{2}-y_{\min}^{2}} = u\,\Delta$ identity. The shipped utility uses the *quadratic-in-* $y^{2}$ variant derived above, which is globally exact and unifies the two endpoint patterns under a single weight formula. Phase 5 itself was retreated as research-grade (continuous-µ $K_{\rm bc}$ is hypersingular and cannot be wired into production); the substitution is the durable promotion-worthy primitive recovered from the wreckage. See §22.5 for the Gauss-Jacobi fallback when a *both-endpoint* $\sqrt{}$-vanishing integrand makes the visibility-cone substitution insufficient, and `chord_half_lengths` for the chord primitive whose annular partition the substitution is built to complement.

</details>

---

## Investigator memos

- `.claude/agent-memory/numerics-investigator/phase5_*.md` (5 documents covering each round).
- Literature memo at `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md`.
- Original diagnostic at `derivations/diagnostics/diag_phase5_round3_visibility_cone_quad.py` (uses the original linear-shift variant; the shipped utility uses the global quadratic-in-$y^2$ variant per §22.7).

References:
- [SanchezTTSP1986] R. Sanchez (1986), *Transport Theory & Statistical Physics* 14, 333–343, DOI 10.1080/00411458608210456.

Closing posture: Issue #133 closed wontfix. The continuous-µ $K_{\rm bc}$ cannot be wired into production (hypersingular kernel; no Nyström-compatible discretisation exists). For multi-bounce specular today, use `boundary="specular_multibounce"` (Phase 4 matrix-Galerkin form). The visibility-cone substitution is the durable promotion-worthy primitive recovered from the investigation; it lives in the production `chord_quadrature` recipe per §22.9. If a future session resurrects continuous-µ work (e.g. via a different singular-integral-equation framework — Hadamard finite-part with proper regularisation, or a fundamentally different operator), a new tracking issue is needed.
