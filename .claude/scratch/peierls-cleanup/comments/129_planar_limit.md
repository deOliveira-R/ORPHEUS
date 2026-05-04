<!-- Posted as part of the post-#138 Peierls documentation cleanup, 2026-04-30 -->
<!-- Source: docs/theory/peierls_nystrom.rst lines 885-907 (~22 LoC) -->
<!-- Cleanup commit: <!-- COMMIT-HASH --> -->

## Issue #129 investigation log — Planar-limit cross-check (hollow cylinder → slab)

**Summary (2026-04-30 — relocated from docs to issue):**

1. The Phase G.4 plan originally hoped that **a hollow cylinder at $r_0 \to R$ would reduce to a slab of thickness $L = R - r_0$ at the $10^{-8}$ level** as a cross-geometry verification cross-check. **It does not.**
2. Empirical probe at $r_0 = 0.999\,R$, $R = 1$, $L = 10^{-3}$, $\Sigma_t = 1$, $c = 0.4$, $\nu\Sigma_f = 0.6$: unified slab gives $k_{\rm eff} = 0.002\,355$, unified cylinder gives $k_{\rm eff} = 0.001\,825$ — a **22 % relative disagreement** at matched optical parameters.
3. The physical reason is **not a bug**: the cylinder's $\mathrm{Ki}_1$ has already integrated the axial direction analytically, giving a different chord-length distribution than the slab's $E_1$. There is no simple thin-shell equivalence at matched $(\Sigma_t, \Sigma_s, \nu\Sigma_f, L)$. This issue tracks the active research record; previously sat as an OQ bullet inside §slab-polar with no issue number, and is being relocated to keep #129's investigation colocated with the issue.

---

<details><summary>Full investigation record (relocated from docs/theory/peierls_nystrom.rst:885–907 — §slab-polar related-OQs)</summary>

### Empirical probe

The naive claim "hollow cylinder at $r_0 \to R$ reduces to a slab of thickness $L = R - r_0$" does **not** hold at the $10^{-8}$ level that Phase G.4 originally hoped for.

Probed empirically at:

- $r_0 = 0.999\,R$
- $R = 1$
- $L = 0.001$
- $\Sigma_t = 1$
- $c = 0.4$
- $\nu\Sigma_f = 0.6$

Result:

- unified slab: $k_{\rm eff} = 0.002\,355$
- unified cylinder: $k_{\rm eff} = 0.001\,825$
- **22 % relative disagreement** at matched optical parameters.

### Physical explanation: chord-distribution mismatch

The physical reason is that the cylinder's $\mathrm{Ki}_1$ has already integrated the axial direction analytically; the remaining in-plane chord-length distribution in a thin annular shell scales as $\sqrt{2\,R\,L}$ for tangential rays ($\approx 0.045$ for the probe's parameters), not $L/|\mu|$ as in a slab.

The two kernels therefore see **different optical-depth spectra** even in the thin-shell limit, so there is no simple geometric equivalence at matched $(\Sigma_t, \Sigma_s, \nu\Sigma_f, L)$. A meaningful planar limit needs either:

1. **A ray-distribution-matched comparison** — explicitly construct slab and cylinder problems whose chord-length PDFs agree, not whose $L$ values agree, or
2. **A curvature-over-thickness expansion** — an asymptotic expansion of the cylinder kernel in $L/R$, identifying the leading-order term as the slab kernel and accounting for the $O(L/R)$ correction explicitly.

Both are future work for this issue. Phase G.4 as specified in the original plan is filed as this GitHub Issue for future physics investigation rather than a shipping test.

</details>

### Why this is being moved out of the doc

Per the post-#138 cleanup directive: **"failed-experiment narrative belongs in GitHub issues, not in evergreen theory docs."** The 22 % disagreement is not a failed experiment in the bug sense — it is correct physics that falsified the naive Phase G.4 plan. The active research record (the chord-distribution explanation, the empirical numbers, and the two candidate paths forward) belongs here in #129; the Sphinx OQ bullet at the slab-polar §related-OQs anchor now reads:

> **OQ — Planar-limit cross-check against hollow cylinder?** Probed empirically and found structurally non-trivial; tracked in [Issue #129](https://github.com/deOliveira-R/ORPHEUS/issues/129). The cylinder's $\mathrm{Ki}_1$ has already integrated axial direction, giving a different chord-length distribution than slab — there is no simple thin-shell equivalence at matched optical parameters.

This way a future agent reading the theory page sees the warning and the pointer here for the empirical evidence + candidate paths forward.

### Cross-links

- Surviving Sphinx OQ: `theory-peierls-slab-polar` §"Related open questions" (one-bullet pointer).
- Cleanup commit: <!-- COMMIT-HASH -->
- Phase G.4 origin: planar-limit cross-check originally specified as part of slab-polar verification rollout.
