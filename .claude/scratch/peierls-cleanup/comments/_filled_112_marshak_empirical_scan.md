# Issue #112 investigation log — empirical canonical-Marshak variant scan (post-#138 docs cleanup)

**Status:** OPEN — Phases A (sphere plateau) and C (cylinder divergence) of the rank-N BC fix remain deferred. A 2026-04-18 empirical scan tested four candidate "canonical" variants against the shipped V1 implementation and found that none of the natural elegant rewrites improve convergence uniformly. The scan results are load-bearing evidence for the active investigation and live here as the durable research record.
**Date:** 2026-04-30
**Why this comment exists:** Per the post-#138 docs cleanup directive, failed-experiment narrative is being relocated from `docs/theory/peierls_nystrom.rst` into the GitHub issue that originated each investigation. This comment is the **active investigation log** for the canonical-Marshak variant scan: the V1/V2/V4/V6 sphere and cylinder $k_{\rm eff}$ error tables (`peierls_nystrom.rst:11117–11205`) plus observations 1–3 (`peierls_nystrom.rst:11209–11223`). The Sphinx page now keeps observation 4 only as a single sentence inside the surviving "What this teaches" subsection; observations 1–3 + the variant-scan tables are the load-bearing empirical evidence for #112 and live here.

[`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)

## Empirical variant scan (verbatim)

Bare homogeneous 1G 1-region white-BC eigenvalue ($\Sigma_t = 1$, $\Sigma_s = 0.5$, $\nu\Sigma_f = 0.75$, $k_\infty = 1.5$):

### Sphere $k_{\rm eff}$ error vs rank $N$

| Variant | Description | R=1 N=1 | R=1 N=2 | R=1 N=3 | R=1 N=8 | R=10 N=1 | R=10 N=8 |
|---|---|---|---|---|---|---|---|
| V1 | Current (shipped) | 27 % | **1.2 %** | 2.3 % | 2.5 % | 0.28 % | 0.17 % |
| V2 | Jacobian mode 0, $R=\mathrm{diag}(2n+1)$ | 50 % | 29 % | 25 % | 24 % | 5.3 % | 5.2 % |
| V4 | Jacobian mode 0, $R=B^{-1}$ | 29 % | 16 % | 15 % | 15 % | 5.2 % | 5.2 % |
| V6 | Cosine-$\mu$ integrand, $R=B^{-1}$ | **1.1 %** | 19 % | 19 % | 19 % | 6.3 % | 6.9 % |

### Cylinder $k_{\rm eff}$ error vs rank $N$

| Variant | Description | R=1 N=1 | R=1 N=2 | R=1 N=3 | R=1 N=8 | R=10 N=1 | R=10 N=8 |
|---|---|---|---|---|---|---|---|
| V1 | Current (shipped) | 21 % | 8.3 % | 27 % | 107 % | 1.1 % | 0.9 % |
| V2 | Jacobian mode 0, $R=\mathrm{diag}(2n+1)$ | 41 % | 17 % | **0.45 %** | 65 % | — | — |
| V4 | Jacobian mode 0, $R=B^{-1}$ | 23 % | 2.6 % | 13 % | 78 % | — | — |

(The V6 cylinder scan was not run; only V1/V2/V4 cylinder data is available.)

## Observations 1–3 (verbatim)

1. **V1 is empirically the best overall.** Sphere converges to a $\sim 2.5\,\%$ plateau at thin $R$; thick cells are sub-percent. Cylinder rank-1 matches legacy Mark accurately.

2. **V2 cylinder at rank-3 R=1 MFP hits the canonical DP_2 prediction** ($0.45\,\%$). This *proves* that the factored structure $G R P$ is capable of canonical rank-$N$ convergence — but only at the cost of a **degraded rank-1** convention (41 % error at rank 1).

3. **V6's cosine-integrand rank-1 is accidentally good for sphere at R=1** (1.1 %) but its rank-$N$ never improves. The rank-1 match is a geometric coincidence at thin $R$ where $\mu_{\rm exit} \approx 1$ for most rays, not a canonical result.

## What the scan teaches (active investigation context)

The investigation reveals a **conceptual entanglement** between the mode-0 convention and the rank-$N$ structure: the legacy mode-0 convention (isotropic-source escape probability, no Jacobian, paired with $R_{00} = 1$) is *not* the canonical partial-current moment, and switching to a canonical mode-0 convention (cosine-weighted moment with $(\rho_{\max}/R)^2$ Jacobian or explicit $\mu_{\rm exit}$ weight) gives a *different* rank-1 result. None of the four candidate rewrites tested in this scan (V2 / V4 / V6, plus the implicit V1 baseline) decouples the two, hence the deferred status — the surviving Sphinx narrative in §29 of `peierls_nystrom.rst` retains the continuous-Marshak derivation, the conceptual-entanglement lesson, and the Stepanek path-forward but no longer carries the variant-scan tables (those are this comment).

## Diagnostics provenance

The scan was performed on 2026-04-18 via the diagnostics scripts:

- `diag_rank_n_09_*.py` through `diag_rank_n_12_*.py` — the V1/V2/V4/V6 sphere + cylinder runs.
- `diag_rank_n_13_phaseAC_summary.md` — the prose summary that was lifted into §29 of `peierls_nystrom.rst` (and which this comment now carries the load-bearing tables out of).

The scripts live under `derivations/diagnostics/` and remain as reproducible artifacts for the next Phase A / Phase C session.

## Cross-links

- Sphinx (post-cleanup): `docs/theory/peierls_nystrom.rst` §29 "Phases A and C: canonical Marshak investigated, deferred" retains intro + continuous-Marshak derivation + the conceptual-entanglement lesson + Stepanek path-forward + deferred-work status block. Observation 4 ("Changing $R$ from $\mathrm{diag}(2n+1)$ to $B^{-1}$ does not rescue a degraded mode-0 convention") survives as a one-sentence summary.
- Code: `orpheus.derivations.peierls_geometry.compute_P_esc_mode_marshak` (the dead-code Phase A/C entry point retained for the future fix; see also Issue #119/#122/#123 sibling rank-N infrastructure threads).
- Issues: this issue (#112) is the canonical home for the Phase A/C deferred work; #122 (Lambert/Marshak gauge) carries the gauge-ambiguity investigation that motivated the variant scan; #123 (rank-N stability protocol L19) carries the stability infrastructure.
- Three candidate rewrites that the scan tested (verbatim from §29 continuous derivation):
  1. **Mode-0 convention**: use the canonical cosine-weighted moment for mode 0 as well as for $n \ge 1$ (instead of the legacy isotropic-source escape probability).
  2. **Reflection operator**: replace $R = \mathrm{diag}(2n+1)$ with $R = B^{-1}$.
  3. **Integrand form**: use $\mu_{\rm exit}$ as the explicit cosine weight in the integrand, rather than the $(\rho_{\max}/R)^2$ surface-to-observer Jacobian.

These map to the V2 / V4 / V6 columns in the scan tables above.
