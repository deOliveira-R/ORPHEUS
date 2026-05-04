# Issue #122 close-out — Lambert/Marshak gauge investigation (Direction Q)

**Date:** 2026-04-30
**Doc cleanup commit:** [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
**Source relocated from:** `docs/theory/peierls_nystrom.rst` lines 4580–4793.

---

## Summary

- F.4 scalar (Hébert 2009 Eq. 3.323) pairs **Lambert-basis** $P_{\rm esc}, G_{\rm bc}$ with **Marshak-basis** $W$. The pairing is formally inconsistent (different inner products) yet empirically load-bearing: every formally-consistent rank-N closure we built plateaus 10–100× above F.4's residual on the hollow sphere at $\sigma_t R \ge 5$.
- **Verdict (B): F.4's Lambert/Marshak asymmetry is a lucky rank-1 algebraic accident with a principled re-phrasing.** At rank-1 the mismatch is a scalar gauge ($M^{(1)} = \sqrt{2}/2$) that factors out of the closure; at rank $N \ge 2$ the basis change is a *genuine basis rotation* (singular values of $M^{(2)}$ are unequal, $(0.460, 0.888)$) that no scaling vector $\beta_{\rm eff}$ can absorb. The off-diagonal entries of $M$ propagate through $(I - \beta W)^{-1}$ with uncontrolled amplification — this is the algebraic origin of the rank-N Lambert-P/G catastrophe (33–737 % $k_{\rm eff}$ error on the 6-point reference grid).
- **Production decision (canonical, preserved): the guard on `boundary="white", n_bc_modes > 1` stays in place. F.4 remains production.** The rank-N Marshak primitives are retained because they are tested, reciprocity-verified, and equal-to-F.4 at $N = 1$ — but no principled rank-N generalisation of the Lambert-side trick exists. Documented in the surviving §F.5 production-decision stub (label `peierls-rank-n-per-face-closeout`).

Surviving doc anchors:
- `:label:` `peierls-half-range-inner-products` — definition of the L and M inner products.
- `:label:` `peierls-change-of-basis` — definition of $M_{nm}$ and the Lambert→Marshak expansion.
- `:label:` `peierls-M-rank-1` (scalar gauge $\sqrt{2}/2$).
- `:label:` `peierls-M-rank-2` (upper-bidiagonal genuine rotation).
- `:label:` `peierls-WM-WL-asymmetric` (the asymmetric $W_M = B^\mu W_L$ identity at infinite rank).
- `:ref:` `peierls-f4-rank-1-gauge-why` (the section title that survives as a stub pointing here).
- `:ref:` `peierls-rank-n-stability` (Issue #123 — the L19 protocol any future candidate must clear).

<details><summary>Full investigation record (verbatim from peierls_nystrom.rst:4580–4793)</summary>

### Why F.4 works at rank-1 but does not generalise (Direction Q, Issue #122)

F.4 scalar (Eq. 3.323 of Hébert 2009) pairs **Lambert-basis** $P_{\rm esc}, G_{bc}$ (integrand has no $\mu$ weight on the outgoing side) with **Marshak-basis** $W$ (the transmission operator is defined on the $\mu$-weighted half-range Gram). The pairing is formally inconsistent — the escape / coupling primitives and the transmission operator live in different inner products — yet it is **empirically load-bearing**: every formally-consistent rank-N closure we built (Marshak everywhere, Lambert everywhere, split $c_{\rm in}$-aware basis with adaptive scale — see `diag_cin_aware_split_basis_keff.py`) plateaus 10–100× above F.4's residual on the hollow sphere at $\sigma_t R \ge 5$.

Two sub-agent investigations closed Issue #122 with verdict **(B)**: **F.4's Lambert/Marshak asymmetry is a lucky rank-1 algebraic accident with a principled re-phrasing.** The algebraic bridge between F.4 and the rank-N Marshak closure is the solid-harmonic change-of-basis matrix between the two half-range $L^2$ structures on $[0, 1]$.

### Setup

On the outgoing-$\mu$ half-line the two natural inner products are

```math
\langle f, g \rangle_L \;=\; \int_0^1 f(\mu)\,g(\mu)\,\mathrm d\mu,
\qquad
\langle f, g \rangle_M \;=\; \int_0^1 f(\mu)\,g(\mu)\,\mu\,\mathrm d\mu.
```
(label: `peierls-half-range-inner-products`)

Shifted-Legendre on $[0, 1]$ is orthogonal under $\langle \cdot, \cdot \rangle_L$; F.4 uses the $L$-orthonormal basis $\phi_n^L(\mu) = \sqrt{2n+1}\,P_n(2\mu-1)$ for its $P$ and $G$. The ORPHEUS rank-N Marshak closure (helper `_build_closure_operator_rank_n_white`, guarded behind `NotImplementedError`) uses the $M$-orthonormal basis $\psi_n^M$ obtained by Gram-Schmidt of $\{1, \mu, \mu^2, \ldots\}$ under $\langle \cdot, \cdot \rangle_M$ (equivalent to scaled Jacobi $P_n^{(0,1)}(2\mu-1)$).

### Change-of-basis matrix

The matrix that expresses Lambert ONB functions in Marshak ONB coordinates is

```math
M_{nm} \;=\; \langle \psi_n^M, \phi_m^L \rangle_M
        \;=\; \int_0^1 \psi_n^M(\mu)\,\phi_m^L(\mu)\,\mu\,\mathrm d\mu,
\qquad
\phi_m^L \;=\; \sum_n M_{nm}\,\psi_n^M.
```
(label: `peierls-change-of-basis`)

The Marshak Gram matrix of the Lambert basis (already cited above for the Phase F.5 rank-N primitives), $B^{\mu}_{mn} = \langle \phi_m^L, \phi_n^L \rangle_M$, equals $M^{\!\top} M$ identically. This is the $B^{\mu}$ matrix consumed by the existing `_build_closure_operator_rank_n_white` helper.

### Symbolic closed forms

Verified in `derivations/diagnostics/diag_lambert_marshak_basis_change.py` (runs in ~1 s).

At rank $N = 1$ the change-of-basis matrix is a **scalar**

```math
M^{(1)} \;=\; \tfrac{\sqrt{2}}{2} \;\approx\; 0.7071,
\qquad
(B^{\mu})^{(1)} \;=\; \tfrac{1}{2}.
```
(label: `peierls-M-rank-1`)

At rank $N = 2$ the change-of-basis matrix is **upper bidiagonal**

```math
M^{(2)} \;=\;
\begin{pmatrix}
  \tfrac{\sqrt{2}}{2} & \tfrac{\sqrt{6}}{6} \\
  0                   & \tfrac{\sqrt{3}}{3}
\end{pmatrix},
\qquad
(B^{\mu})^{(2)} \;=\;
\begin{pmatrix}
  \tfrac{1}{2} & \tfrac{\sqrt{3}}{6} \\
  \tfrac{\sqrt{3}}{6} & \tfrac{1}{2}
\end{pmatrix}.
```
(label: `peierls-M-rank-2`)

$M^{(2)}$ has singular values $(0.460, 0.888)$. They are **unequal**, so $M^{(2)}$ is not a scalar multiple of an orthogonal matrix — it is a **genuine basis rotation**, not a scalar gauge. The rank-3 and rank-4 change-of-basis matrices preserve this upper-bidiagonal structure; the sum of magnitudes of strictly-off-diagonal entries grows monotonically (0, 0.408, 0.855, 1.318 at ranks 1–4).

### Implication for F.4

At rank-1, the Lambert / Marshak mismatch is the scalar $M^{(1)} = \sqrt{2}/2$; it factors out of the closure $K_{bc} = G \beta (I - \beta W)^{-1} P$ as a rescaling of $\beta$. F.4's "trick" is to effectively use the Marshak-side $\beta_{\rm eff} = \beta / (M^{(1)})^2 = 2\beta$ without explicitly saying so. The closure-level gauge factor $\alpha(\tau) \approx 0.38$ measured in `diag_lambert_marshak_symbolic.py` (see §Experiment 7 of the research log) is this basis-change scalar **times** the $\exp(-\sigma R)$ attenuation integrated through the ray — a single number, identifiable, and absorbable.

At rank $N \ge 2$, $M$ ceases to commute with any diagonal $\operatorname{diag}(\beta_0, \beta_1, \ldots)$. There is **no vector** $\beta_{\rm eff}$ that, pre-multiplied onto the Marshak closure, reproduces F.4's Lambert behaviour. The off-diagonal entries of $M$ mix Lambert mode-0 into Marshak mode-1 and vice versa, and that mixing propagates through the closure operator $(I - \beta W)^{-1}$ with uncontrolled amplification. This is the algebraic origin of Experiment E2.4's rank-N Lambert-P/G catastrophe (33–737 % $k_{\rm eff}$ error on the 6-point reference grid) documented in the research log.

### The precise obstruction: asymmetric µ-multiplication with a polynomial-truncation leak

A frame-covariant rewrite of F.4 — conjugating the Lambert closure by $M$ on both sides as

$$K_{bc}^{F.4} = G_L \cdot M^{\!\top}\,(I - M W_L M^{\!\top})^{-1}\,M\,P_L$$

— was tested at rank-2 and rank-3 (symbolic + mpmath, see `derivations/diagnostics/diag_frame_4_connection_form.py`) and fails: 22 % relative error at the rank-2 anchor, 48 % at rank-1, across 5 conjugation variants and 7 values of $\sigma_t R$. $W$ is NOT a (1,1) tensor under $M$; at rank-1 already $M W_L M^{\!\top} = 0.05$ versus $W_M = 0.005$ at $\sigma_t R = 10$, a tenfold discrepancy. The correct algebraic relationship is the **asymmetric identity**

```math
W_M \;=\; B^{\mu}\,W_L,
\qquad
B^{\mu} \;=\; M^{\!\top} M,
```
(label: `peierls-WM-WL-asymmetric`)

verified bit-exactly at infinite rank. At rank $N$, this identity holds on rows $0, 1, \ldots, N-2$ of $W_M - B^{\mu} W_L$ (vanishing symbolically) but **row $N-1$ carries a non-vanishing τ-dependent polynomial-truncation residual** — because $\mu \cdot \tilde P_{N-1}$ has a $\tilde P_N$ component that the rank-$N$ basis cannot represent. The above equation is the exact statement of Lambert ↔ Marshak at the transmission operator; the rank-1 scalar-gauge picture is the finite-truncation limit where the leaking row has empty support. **This obstruction is structural and cannot be cured by any basis rotation on the rank-N space.**

### Connection to the gauge-theoretic literature

Sanchez (2014) NSE 177(1), DOI `10.13182/NSE12-95`, establishes that the first-order $P_N$ equations are **degenerate** — multiple IC/BC closures are admissible, with uniqueness imposed via second-order-parity equivalence and solid-harmonic expansions that reproduce results originally due to Davison (1957) and Rumyantsev (1950s). The Sanchez theorem applies to the differential $P_N$ equations, not directly to the integral CP operator F.4 lives on, but the gauge-freedom framing is the right language: **at rank-1 the integral CP admits a scalar gauge (eq. `peierls-M-rank-1`) that F.4 exploits**; at rank-N the "gauge" becomes a full basis rotation (eq. `peierls-change-of-basis`) that the Marshak closure must account for explicitly. Canonical open-access channels to the underlying solid-harmonic material are Davison (1957) *Neutron Transport Theory* (Oxford; AEC NAA-SR-3509) and Case & Zweifel (1967) *Linear Transport Theory*.

### Statistical-mechanical picture (partial)

The scalar gauge $\alpha(\tau, \rho) \approx 0.38$ has a partial statistical-mechanical interpretation via the surface Markov chain on the outer hemisphere: state space is outgoing $\mu \in [0, 1]$; the transition kernel is ballistic chord through the hollow sphere followed by isotropic re-emission with in-shell scattering $c = 1/3$. Direct 500 k-history Monte Carlo of this chain (see `derivations/diagnostics/diag_frame_3_surface_markov_mc.py`) gives a Perron eigenfunction $p_{\infty}(\mu)$ whose moments are **ρ-independent to $\le 1.3\,\%$** across $\rho \in [0.3, 0.5]$ at fixed $\tau$ — this is the mechanism behind the observed ρ-flatness of $\alpha$. But $p_{\infty}(\mu)$ is **not Laplace-type**: an exponential fit $A e^{-\lambda\mu} + B$ leaves 7–11 % residual; the histogram is monotone *increasing* in $\mu$, with $\mathbb E[\mu] \approx 0.56\text{–}0.61$ and $\mathbb E[\mu^2]/\mathbb E[\mu] \approx 0.70\text{–}0.72$. No natural moment of $p_{\infty}$ identifies $\alpha$ to better than 5 % uniformly across the 6-point grid. The rank-N polynomial expansion of $p_{\infty}$ is therefore **basis-resistant** because $p_{\infty}$ is neither polynomial nor single-exponential — supplementing the algebraic Schur-reduction story (eq. `peierls-WM-WL-asymmetric`) with an independent statistical obstruction. An analytical computation of $p_{\infty}$ as the left Peierls-kernel eigenvector on the outer surface would settle the quantitative $\alpha$ identification without MC bias; this is unresolved and listed as a follow-up.

### Production decision

The guard on `boundary="white", n_bc_modes > 1` stays in place. F.4 remains production. The rank-N Marshak primitives below are retained because they are tested, reciprocity-verified, and equal-to-F.4 at $N=1$ — but no principled rank-N generalisation of the Lambert-side trick exists. Any future rank-N white-BC closure candidate must compete against F.4 under the two-quadrature stability protocol (Direction N, Issue #123) — see L19 in the research log.

</details>

---

## Empirical scan and Stepanek 1981 calibration plan (research record)

The empirical $\alpha(\tau)$ scan and the Stepanek 1981 calibration plan referenced in the source text (Experiment 7 of `.claude/plans/rank-n-closure-research-log.md`) probed the gauge identifiability question. The scan found $\alpha \approx 0.38$ ρ-independent at fixed $\tau$ (within 1.3 %) across the 6-point reference grid, but no closed-form analytical channel back to a half-range polynomial moment of the Perron eigenfunction was found. Stepanek 1981 was identified as a candidate calibration source for the surface partial-current redistribution; the canonical decision archived in the surviving §F.5 production-decision stub is to **not** chase a Stepanek-calibrated rank-N gauge, on the strength of the eq. `peierls-WM-WL-asymmetric` row-$N-1$ truncation residual (the obstruction is structural, not calibration-dependent).

References:
- Hébert (2009/2020), *Applied Reactor Physics* (3rd ed.), Eq. 3.323.
- Sanchez (2014), *Nuclear Science and Engineering* 177(1), DOI `10.13182/NSE12-95`.
- Davison (1957), *Neutron Transport Theory* (Oxford; AEC NAA-SR-3509).
- Case & Zweifel (1967), *Linear Transport Theory*.
- Stepanek (1981) — referenced in the research log as the candidate gauge calibration source.

Closing posture: Issue #122 closed with verdict (B); the Lambert/Marshak gauge is identified, characterised, and shown to be unabsorbable at rank $N \ge 2$. F.4 stays production behind the `n_bc_modes > 1` guard. Future work directs to Issue #123 (L19 stability protocol) and Issue #128 (optional QMC migration of F.4 production quadrature).
