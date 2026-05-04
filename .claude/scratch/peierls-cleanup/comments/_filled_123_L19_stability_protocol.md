# Issue #123 investigation log — Rank-N stability protocol (L19)

**Date:** 2026-04-30
**Doc cleanup commit:** [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
**Source relocated from:** `docs/theory/peierls_nystrom.rst` lines 4794–4889.

**Status:** OPEN. This is the **active research record** — not a post-mortem. The L19 protocol below is the operational gate any future rank-N white-BC closure candidate must clear before it can claim to beat F.4. The next agent picking up rank-N work should treat this comment as the load-bearing definition.

---

## Summary

- A rank-N white-BC closure candidate $C$ claims to beat F.4 at a reference point $(\tau, \rho)$ **only if** the claim survives a two-quadrature **signed-error stability protocol** (S1–S5 below). Single-quadrature comparison rewards quadrature-noise cancellation rather than truncation-residual reduction; this protocol is the operational response to lessons L17–L19 in `.claude/plans/rank-n-closure-research-log.md`.
- Implementation lives at `tests.cp.test_peierls_rank_n_protocol.assert_rank_n_structural_win` and ships with pinning tests at the six reference points $(\sigma_t R, \rho) \in \{5, 10, 20\} \times \{0.3, 0.5\}$. Two of the six (`σ_t R = 10, ρ = 0.3` and `σ_t R = 20, ρ = 0.5`) reproduce the canonical L17 sign-flip of F.4 itself and are tagged `@pytest.mark.slow`.
- **Open work remaining**: the ULTRA = (5, 10, 96) baseline + RICH+pp at all six points is unresolved on devcontainer hardware (exceeded 120-s-per-point budget). Resolving the full ULTRA baseline would make $\mathcal Q$ trivially L19-compliant without needing S5; target is ≥ 300 s per point at $\sigma_t R = 20$ (richer hardware or a relaxed wall budget). See L20 in the research log.
- A **randomised-QMC alternative** (Owen-scrambled Sobol' on the angular dimension, validated 2026-04-22) gives 95 % bootstrap CI widths of $6 \times 10^{-5}$ to $6 \times 10^{-4}$ on F.4 at all six reference points — 20–100× tighter than the PG RICH vs RICH+panels spread the L19 protocol uses to *detect* instability. Both L17 sign-flip points resolve to crisp negative QMC means whose CIs do not cross zero. A future candidate can therefore replace the S3–S4 gates with a **single CI-separation assertion**. Sketched as `assert_rank_n_qmc_structural_win` in the Frame 5 memo; not shipped because no current closure passes Frame 5. Issue #128 tracks the optional migration of F.4 production quadrature from product-Gauss to randomised QMC; LOW priority, not on the critical path.

Surviving doc anchors:
- `:label:` `peierls-rank-n-stability` — the equation block defining S1–S5 (this label survives in `peierls_nystrom.rst`).
- `:ref:` `peierls-rank-n-per-face-closeout` — the surviving §F.5 production-decision stub (F.4 remains production).
- `:ref:` `peierls-f4-rank-1-gauge-why` — Issue #122 sister investigation explaining *why* F.4 is the production reference.

<details><summary>Full L19 protocol (verbatim from peierls_nystrom.rst:4794–4889)</summary>

### Rank-N stability protocol (L19)

A rank-N white-BC closure candidate $C$ claims to beat F.4 at a reference point $(\tau, \rho)$ if and only if the claim survives the **two-quadrature signed-error stability protocol** defined here. This protocol is the operational response to lessons L17–L19 in `.claude/plans/rank-n-closure-research-log.md`: RICH = (4, 8, 64) is below F.4's own structural floor at $\sigma_t R \ge 10$, so a naive single-quadrature comparison rewards quadrature-noise cancellation rather than truncation-residual reduction.

### Protocol

Let $\mathcal Q = (q_1, q_2, \ldots, q_K)$, $K \ge 2$, be a sequence of quadrature triples $(n_{\rm panels}, p_{\rm order}, n_{\rm ang})$ of monotonically increasing refinement (lexicographic in any component). Let $e_C(q_k)$ and $e_{F.4}(q_k)$ denote the **signed** relative errors $(k_{\rm eff}^{C/F.4}(q_k) - k_\infty)/k_\infty$ at each quadrature.

The candidate $C$ is a **structural win** over F.4 at $(\tau, \rho)$ iff all five of the following hold:

```math
\begin{aligned}
& \textbf{(S1)} \quad K \ge 2, \\
& \textbf{(S2)} \quad |e_C(q_k)| \;<\; |e_{F.4}(q_k)|
    \quad \text{for every } k, \\
& \textbf{(S3)} \quad \operatorname{sign}\,e_C(q_k) \;=\; \operatorname{sign}\,e_C(q_{k+1})
    \quad \text{for every } k, \\
& \textbf{(S4)} \quad |e_C(q_k)| \;\ge\; |e_C(q_{k+1})|
    \quad \text{for every } k, \\
& \textbf{(S5)} \quad \operatorname{sign}\,e_{F.4}(q_k)
    \;=\; \operatorname{sign}\,e_{F.4}(q_{k+1})
  \; \wedge \; |e_{F.4}(q_k)| \;\ge\; |e_{F.4}(q_{k+1})|
    \quad \text{for every } k.
\end{aligned}
```
(label: `peierls-rank-n-stability`)

Assertion **(S5)** is the L17/L19 reference-verifiability gate: if F.4 itself sign-flips or grows in magnitude under refinement at $(\tau, \rho)$, the *reference* is unverifiable at $\mathcal Q$ and no rank-N comparison there is admissible. This is strictly stronger than L16's "match quadrature across compared closures" — it additionally demands that the match resolves the smaller structural floor.

### Implementation

The helper `tests.cp.test_peierls_rank_n_protocol.assert_rank_n_structural_win` raises `AssertionError` on any of S1–S5 failing and returns a `StabilityReport` dataclass with the full signed-error trajectory otherwise. The helper ships with pinning tests for the RICH vs RICH+panels pair at the six reference points $(\sigma_t R, \rho) \in \{5, 10, 20\} \times \{0.3, 0.5\}$; two of those six ($\sigma_t R = 10, \rho = 0.3$ and $\sigma_t R = 20, \rho = 0.5$) reproduce the canonical L17 sign-flip of F.4 itself and are tagged `@pytest.mark.slow`.

A baseline of F.4 at ULTRA = (5, 10, 96) and richer at all six points — the reference that would make $\mathcal Q$ trivially L19-compliant without needing S5 — is currently **unresolved** on the devcontainer hardware used during Issue #123 development: ULTRA and RICH+pp exceeded the 120-s-per-point budget at every point tested. Resolving the full ULTRA baseline requires either richer hardware or a relaxed wall budget (target: $\ge$ 300 s per point at $\sigma_t R = 20$). See L20 in the research log.

### Randomized QMC alternative (validated 2026-04-22)

The product-Gauss bias that produces L17 sign flips is driven by the tangent-angle kink in the $\exp(-\tau d)$ integrand, whose Hardy-Krause variation is *bounded in* $\tau$. Owen-scrambled Sobol' on the angular dimension (32 scrambles × 4096 points) gives 95 % bootstrap CI widths of $6 \times 10^{-5}$ to $6 \times 10^{-4}$ per cent on F.4 at all six reference points — 20–100× tighter than the PG RICH vs RICH+panels spread the L19 protocol uses to *detect* instability. Both L17 sign-flip points ($\sigma_t R = 10, \rho = 0.3$ and $\sigma_t R = 20, \rho = 0.5$) resolve to crisp negative QMC means whose CIs do not cross zero. See `derivations/diagnostics/diag_f4_qmc_quadrature.py` and the Frame 5 memo. A future rank-N closure candidate can therefore replace the S3–S4 gates of `peierls-rank-n-stability` by a **single CI-separation assertion**:

> closure mean, F.4 mean with disjoint 95 % CIs, AND closure CI strictly tighter than $|\text{F.4 mean}|$.

The thin wrapper `assert_rank_n_qmc_structural_win(closure, f4, point, N=4096, n_scrambles=32)` implementing this is sketched in the Frame 5 memo; not shipped because no current closure passes Frame 5. Issue #128 tracks the optional migration of F.4 production quadrature from product-Gauss to randomized QMC; LOW priority, not on the critical path.

</details>

---

## Hand-off for the next rank-N attempt

If you are a future agent picking up rank-N white-BC work:

1. **Read this comment in full first.** The L19 protocol is the operational definition of "beats F.4" — any candidate that does not clear S1–S5 (or the Frame 5 QMC CI-separation alternative) is not admissible as production.
2. **Read the sister Issue #122 close-out** for the *why* — the algebraic obstruction (Lambert/Marshak basis rotation at $N \ge 2$) is what every formally-consistent rank-N attempt has hit. Any candidate must either (a) explicitly handle the off-diagonal $M^{(N)}$ coupling, or (b) abandon the Marshak basis and find a different starting structure.
3. **Read Issue #132** (Class B MR×MG falsification) for the heterogeneous-cell version of the same obstruction — the rank-N path was independently falsified there at MR×MG with a ~+57 % sphere catastrophe. The Hébert (1−P_ss)⁻¹ closure is the active partial-fix for sphere only; cylinder and slab still raise NotImplementedError.
4. **Tooling already in place:**
   - `tests.cp.test_peierls_rank_n_protocol` — protocol harness with the six reference-point fixtures.
   - `derivations/diagnostics/diag_f4_qmc_quadrature.py` — QMC reference for Frame 5.
   - `.claude/plans/rank-n-closure-research-log.md` — the full L0–L21 lessons trail.

References:
- The L19 protocol's own definition above.
- L17/L19/L20 in the rank-n-closure research log.
- Frame 5 memo (referenced as `.claude/agent-memory/numerics-investigator/frame_5_qmc_quadrature.md`).
- Issue #128 (optional F.4 → QMC migration; LOW priority).

Closing posture: Issue #123 stays OPEN as the active research record. The protocol itself is operational; what remains open is (a) the ULTRA baseline at all six points, and (b) any candidate closure that actually clears S1–S5 (none does today).
