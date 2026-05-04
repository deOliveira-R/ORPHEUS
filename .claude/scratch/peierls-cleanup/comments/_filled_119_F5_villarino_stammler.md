# Issue #119 close-out — F.5 rank-N per-face: Villarino-Stamm'ler per-mode extension falsified

**Date:** 2026-04-30
**Doc cleanup commit:** [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
**Source relocated from:** `docs/theory/peierls_nystrom.rst` lines 4443–4538.

---

## Summary

- A novel per-diagonal-mode extension of Villarino-Stamm'ler (V-S) Gauss-Seidel a-posteriori normalisation (Hébert 2009 Eqs. 3.347–3.352) was tested as a candidate fix for the **rank-N per-face white-BC plateau at 1.42 %** (hollow sphere, $r_0/R = 0.3$, $\sigma_t \cdot R = 5$, $N \ge 2$).
- **Falsified, 2026-04-21**: enforcing per-mode conservation makes k_eff *worse* on the µ-ortho pipeline (1.42 % → 1.87 %) and is essentially a no-op on the shipped Marshak pipeline (10.86 % → 10.83 %). **The plateau is a cross-mode coupling failure, not a conservation-forcing failure** — correcting only the diagonal blocks distorts the balance between diagonal and off-diagonal coupling.
- **Production decision (preserved): F.4 scalar rank-2 per-face is production.** Stamm'ler IV Eq. 34 rank-2 closure is the production decision (already documented in the surviving §F.5 production-decision stub at `docs/theory/peierls_nystrom.rst` ≈ lines 4540–4578, label `peierls-rank-n-per-face-closeout`).

Surviving doc anchors:
- `:label:` `hebert-3-350` (the V-S identity itself, Hébert p. 129).
- `:label:` `peierls-change-of-basis`, `peierls-M-rank-1`, `peierls-M-rank-2`, `peierls-WM-WL-asymmetric` (Lambert/Marshak basis-rotation algebra that explains why per-mode V-S cannot help).
- `:ref:` `c-in-remapping` (the structural obstruction section that this falsification belongs to).
- `:ref:` `peierls-f4-rank-1-gauge-why` (Direction Q / Issue #122, sister investigation).

---

## Production decision (preserved in Sphinx)

F.4 scalar rank-2 per-face is production. Its residual across optical thickness at $r_0/R = 0.3$ (from `derivations/diagnostics/diag_rank_n_closure_characterization.py`):

| $\sigma_t \cdot R$ | $|k_{\rm eff} - k_\infty| / k_\infty$ | Regime                                                |
|--------------------|---------------------------------------|-------------------------------------------------------|
| 1                  | 3.27 %                                | Streaming limit (quadrature-dominated)                |
| 2.5                | 0.55 %                                | Intermediate                                          |
| 5                  | **0.077 %**                           | Near quadrature floor (Mark DP$_0$ truncation ~0.04–0.1 %) |
| 10                 | 0.26 %                                | Forward-peaked flux                                   |
| 20                 | 0.45 %                                | Strongly absorbing                                    |

The 0.077 % minimum at $\sigma_t \cdot R = 5$ coincides with the Mark DP$_0$ truncation floor for this geometry. This is the expected expression of Stacey 2007's warning on p. 329 that $[E_2(\Sigma)]^N \ne E_2(N\Sigma)$ for DP-0 subdivision of forward-peaked boundary fluxes: the error class is geometric, baked into the DP-0 assumption, and cannot be removed without changing the angular representation.

<details><summary>Full investigation record (verbatim from peierls_nystrom.rst:4443–4538)</summary>

### Villarino-Stamm'ler per-mode extension (novel, falsified 2026-04-21)

Hébert 2009 Ch. 3 explicitly warns (p. 129) that the rank-0 DP$_N$ primitives are *not* guaranteed conservative and recommends an a-posteriori **Villarino-Stamm'ler normalisation** (Eqs. 3.347–3.352). Villarino-Stamm'ler (V-S) is a 30-line Gauss-Seidel fixed-point iteration that multiplies the symmetric rank-0 $t$ matrix by an additive symmetric correction to force row conservation while preserving reciprocity.

```math
\hat{t}_{\ell m} \;=\; (z_{\ell} + z_{m})\, t_{\ell m},
\qquad \ell,\, m \;=\; 1, \ldots, \Lambda + I.
```

Reciprocity is preserved by construction because the scalar factor $(z_{\ell} + z_{m})$ is symmetric in $\ell \leftrightarrow m$ and $t_{\ell m}$ is already symmetric (Hébert p. 129).

Hébert defines V-S only on rank-0 primitives. The novel extension tested here is **per-diagonal-mode V-S on the rank-N W**: for each mode $n \in \{0, \ldots, N{-}1\}$, extract the 2×2 diagonal block

```math
W_{\rm sub}[n] \;=\;
\begin{pmatrix}
   W[n,\, n] & W[n,\, N{+}n] \\
   W[N{+}n,\, n] & W[N{+}n,\, N{+}n]
\end{pmatrix},
```

solve the 2-unknown V-S system for $(z_{\rm outer}^{(n)}, z_{\rm inner}^{(n)})$, and apply the additive symmetric correction per mode. The target is the F.4 scalar row sum (mode-0 conserves against F.4's target as a sanity check; $n \ge 1$ is forced onto F.4's mode-0 target, which is the strongest possible demand for per-mode conservation). Off-diagonal cross-mode blocks $n \ne m$ are left untouched.

#### Per-mode V-S residuals on hollow sphere ($r_0/R = 0.3$, $\sigma_t \cdot R = 5$)

| $N$ | µ-ortho raw | µ-ortho + V-S       | Shipped Marshak raw | Shipped Marshak + V-S |
|-----|-------------|---------------------|---------------------|-----------------------|
| 1   | 2.55 %      | 2.17 %              | 13.53 %             | 13.53 %               |
| 2   | **1.42 %**  | **1.87 % WORSE**    | 10.86 %             | 10.83 %               |
| 3   | 1.42 %      | 1.87 % WORSE        | 10.70 %             | 10.66 %               |
| 4   | 1.43 %      | 1.88 % WORSE        | 10.70 %             | 10.66 %               |

Reciprocity is preserved at machine precision ($A_i\, W_{ij} - A_j\, W_{ji} \le 10^{-16}$) for every row of the table, both pre- and post-V-S, confirming that the additive symmetric form does what Hébert Eqs. 3.350 claim (p. 129). The V-S scheme *also* hits its design target: at $\sigma_t = 0$ the per-mode row sums are driven to the F.4 scalar targets for every mode $n$ (e.g., at $N = 3$: outer row sum $2.137 \to 1.910$, inner $0.118 \to 0.090$ for $n = 0$; outer $1.191 \to 1.910$, $0.547 \to 1.910$ for $n = 1, 2$).

**The falsification is unambiguous.** Enforcing per-mode conservation does not close the plateau — on the µ-ortho pipeline it makes k_eff **worse** (1.42 % → 1.87 %); on the shipped Marshak pipeline it is essentially a no-op (10.86 % → 10.83 %). **The plateau is a cross-mode coupling failure, not a conservation-forcing failure.** Correcting the diagonal blocks distorts the balance between diagonal and off-diagonal coupling, shifting the closure further off rather than toward $k_\infty$.

Diagnostic at `derivations/diagnostics/diag_rank_n_villarino_stammler_per_mode.py`; full memo at `.claude/agent-memory/numerics-investigator/peierls_villarino_stammler_per_mode.md`.

</details>

---

## Why this fits the broader rank-N obstruction

The companion §F.4 / Issue #122 close-out establishes that the Lambert-vs-Marshak basis change is a *genuine basis rotation* at $N \ge 2$ (singular values of $M^{(2)}$ are $(0.460, 0.888)$ — unequal, see eq. `peierls-M-rank-2`). At rank-1 the mismatch is a scalar gauge ($M^{(1)} = \sqrt{2}/2$, eq. `peierls-M-rank-1`) that factors out; at rank $\ge 2$ it does not. The per-mode V-S falsification documented here is the **direct empirical confirmation** of that algebraic argument: the off-diagonal coupling that V-S leaves untouched is exactly where the basis-rotation mass lives, so per-mode conservation cannot help the plateau.

References:
- Hébert (2009/2020), *Applied Reactor Physics* (3rd ed.), Ch. 3 §3.8.5, Eqs. 3.347–3.352 (V-S derivation).
- Stamm'ler & Abbate (1983), *Methods of Steady-State Reactor Physics in Nuclear Design* — Stamm'ler IV Eq. 34 rank-2 closure (production reference).
- Sanchez & McCormick (1982), *Nuclear Science and Engineering* 80(4), 481–535, §III.F.1 Eqs. 165–169 (rank-N partial-current-moment basis).
- Stacey (2007), *Nuclear Reactor Physics*, p. 329 (DP$_0$ subdivision warning).

Closing posture: F.5 rank-N per-face is closed-falsified. Stamm'ler IV Eq. 34 rank-2 closure remains production. Future rank-N white-BC candidates must clear the L19 two-quadrature signed-error stability protocol (Issue #123) before they can claim to beat F.4.
