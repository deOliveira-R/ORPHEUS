# Issue #132 investigation log — Class B rank-N MR×MG: probe-cascade root cause + follow-up directions

**Date:** 2026-04-30
**Doc cleanup commit:** [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
**Source relocated from:** `docs/theory/peierls_nystrom.rst` lines 5267–5409 (Probe G LEGACY/CANONICAL cascade) **plus** lines 5807–5970 (Davison + Augmented Nyström + Cyl Hébert follow-up directions).

**Status:** OPEN. This is the **active research record** — not a post-mortem. The Hébert (1−P_ss)⁻¹ closure is shipped as a substantial partial fix for sphere; cylinder + slab still raise NotImplementedError; the +56 % sphere 1G/2R rank-2 catastrophe is pinned by `xfail-strict=True` regression tests.

**Hard evidence pin:** ERR-030 in `tests/l0_error_catalog.md` — the bug catalog entry (failure-mode classification, how it hid, lesson).

---

## Summary

- **The catastrophe.** Sphere 1G/2R at rank-2 produces $(k_{\rm eff} - k_\infty)/k_\infty = +56.7\,\%$ vs $-1.10\,\%$ on the 1R control; the 2R configuration plateaus to ~+67 % at high $N$ (i.e. *increasing* $N$ makes the answer monotonically worse). Cylinder 1G/2R rank-2 shows the analogous +18 % catastrophe. RICH-preset re-check confirms structural (0.022 % shift between BASE and RICH).
- **Root cause (Probe G localisation).** The legacy `compute_P_esc` / `compute_G_bc` at mode 0 omit the $\rho_{\max}^2$ Jacobian that modes $n \ge 1$ carry, putting mode-0 in the **Lambert** inner product and modes $n \ge 1$ in the **Marshak** inner product. At rank-1 the mismatch is the scalar gauge $M^{(1)} = \sqrt{2}/2$ (eq. `peierls-M-rank-1`); at rank $N \ge 2$ it becomes a *genuine basis rotation* (eq. `peierls-M-rank-2`) and the mode-0 / mode-$n \ge 1$ hybrid is structurally inconsistent.
- **Why pure-canonical is not the fix.** Switching mode-0 to `compute_P_esc_mode(n=0)` makes the rank-N expansion internally consistent but converges to the wrong limit (~−25 % across all 1R/2R configurations at high $N$): the canonical Marshak partial-current basis with mode-0 weighted by $\rho_{\max}^2$ does not reduce to the production rank-1 Mark closure at $N = 1$. It is a *different (worse)* single-mode closure, breaks the rank-1 regression gate, and converges wrong.
- **Production decision (preserved):** Class B ships only the rank-1 Mark closure (`solve_peierls_*_1g` with default `n_bc_modes=1`). Calls with `n_bc_modes ≥ 2` are reachable but **not safe**. The XFAIL-strict tests at `tests/derivations/test_peierls_rank_n_class_b_mr_mg.py` pin the catastrophe magnitude; they flip from xfail-pass to unexpected-pass when Issue #132 lands a corrective re-derivation. **The Hébert (1−P_ss)⁻¹ closure** (`boundary="white_hebert"`) ships as the substantial partial fix for sphere — see surviving §`peierls-class-b-sphere-hebert`.

The doc keeps the core empirical evidence at `peierls_nystrom.rst` ≈ lines 5109–5217 (sphere/cyl signed-error tables; three observations) and the production decision at lines 5410–5476. This comment relocates only the cascade-narrative + follow-up-direction prose.

Surviving doc anchors:
- `:label:` `peierls-class-b-Jn-canonical` — the canonical $J^+_n$ formula (eq. derived from `peierls-rank-n-jacobian-derivation`).
- `:label:` `peierls-rank-n-jacobian-derivation` — the rank-N Jacobian derivation that the Probe G analysis cites.
- `:label:` `peierls-cyl-Pss-derivation`, `peierls-cyl-Gbc-3d-final` — cylinder Hébert primitives (the "Cylinder Hébert" follow-up direction is *partially shipped* in code).
- `:ref:` `peierls-class-b-sphere-hebert` — the Hébert (1−P_ss)⁻¹ partial-fix section that survives.
- `:ref:` `peierls-f4-rank-1-gauge-why` — Issue #122 sister investigation (same Lambert/Marshak basis-rotation algebra).

---

## Probe-cascade narrative

<details><summary>Full forensic detail — Probe B/C/D/E/F localisation cascade (verbatim from peierls_nystrom.rst:5267–5409)</summary>

### Root cause — Probe G normalisation mismatch

The full probe-cascade (Probes B–H, see `derivations/diagnostics/diag_class_b_rank_n_probe_{b,c,d,e,f,g}_*.py`) ruled out:

- **Volume kernel multi-region path** (Probe B, `vacuum_2r`) — the vacuum-BC $K_{\rm vol}$ MR routing is tight against the 1R baseline within $\sim 2 \times 10^{-4}$. The volume kernel is not the bug.
- **Routing path under uniform $\Sigma_t$** (Probe C, `homogeneous_2r`, promoted to permanent passing test `test_class_b_mr_routing_invariance_uniform_sigma`) — sphere/cyl with `radii=[0.5, 1.0]` and uniform $\Sigma_t = 1$ matches `radii=[1.0]` within $5 \times 10^{-3}$ (Issue #114 noise floor). The 2R routing path is consistent with 1R at the $\sim 10^{-3}$ level when the $\Sigma_t$ profile is flat — *the divergence requires a real $\Sigma_t$ step*.
- **Primitive convergence under quadrature refinement** (Probe D, `primitive_quadrature`) — `compute_P_esc_mode` / `compute_G_bc_mode` plateau at $\sim 10^{-5}$ under `n_angular` refinement to 192. The primitives are essentially correct; no closed-form-avoidance anti-pattern à la Issue #131.
- **Conservation defect localisation** (Probe E, `conservation`) — per-node $(K \cdot 1 - \Sigma_t)/\Sigma_t$ defects are 5–7 % rms in 2R-Z, but the 1R control has 9 % rms defect with $k_{\rm eff}$ *still right*. **The conservation defect is not a strong predictor of $k_{\rm eff}$ error in MR.** This is important methodologically: the conservation row-sum test at `tests/derivations/test_peierls_rank_n_conservation.py` uses uniform $\Sigma_t = 1$ where $K \cdot \mathbf 1 = \mathbf 1$ holds by construction; **the test is structurally blind to MR mismatches** (the identity becomes an integrated identity, not pointwise, when $\Sigma_t$ is piecewise). The numerics-investigator's initial conservation diagnostic was structurally wrong for this reason; the lesson is captured in ERR-030 and pinned by the new MR routing-invariance test.
- **Per-mode K_bc isolation** (Probe F, `mode_isolation`) — adding mode-1 alone to mode-0 jumps $k_{\rm eff}$ by **+84 %** on sphere 1G/2R (vs +35 % on the 1R control). Mode-1 contribution does not scale linearly between 1R and 2R — the mode-1 → mode-0 ratio depends sensitively on the $\Sigma_t$ profile in a way that no per-mode primitive could be wrong about (it would need to go wrong at *first contact* with mode 1).

Probe G (`normalization_mismatch`) is the localisation. It runs the same 1R / 2R sweep with two variants of the mode-0 routing:

- **LEGACY** (production, the bug): `compute_P_esc` + `compute_G_bc` at mode 0; no $(\rho_{\max}/R)^2$ Jacobian. Modes $n \ge 1$ use `compute_P_esc_mode` + `compute_G_bc_mode` (with Jacobian).
- **CANONICAL**: `compute_P_esc_mode(n=0)` + `compute_G_bc_mode(n=0)` at mode 0, identical Jacobian-weighted form as modes $n \ge 1$.

#### Probe G — sphere $(k_{\rm eff} - k_\infty)/k_\infty$ at rank-2, BASE preset

| Configuration | LEGACY (production, mismatched) | CANONICAL (consistent, all-Jacobian) |
|---|---|---|
| 1R control ($\Sigma_t = 1$, $k_\infty = 1.5$) | **−1.10 %** | −29.3 % |
| 2R-Z ($\Sigma_t = [1, 2]$, $k_\infty = 0.648$) | **+56.7 %** | −28.0 % |

The CANONICAL variant produces *consistent* errors across 1R and 2R (both plateau near −25 % to −29 % at high $N$) — the mismatch is gone, the routing is internally consistent. **But the CANONICAL closure is uniformly worse**. The legacy mode-0 form is not a *bug* in the sense of "wrong code"; it is a *calibration*: the legacy `compute_P_esc` was historically tuned to make the rank-1 Mark closure (no rank-N involved) approximately right on solid cells. When summed with mode-$n \ge 1$ to form a rank-N outer-product expansion, that calibration breaks the algebraic structure of the rank-N partial-current basis because the two terms do not live in the same expansion space.

### Why the legacy/canonical hybrid is structurally inconsistent

The Marshak partial-current moment of order $n$ from a uniform unit volumetric source at radial node $r_i$ is, by `peierls-rank-n-jacobian-derivation`,

```math
J^{+}_n(r_i) \;=\; \frac{1}{A_d}\,
                   \int_\Omega \tilde P_n(\mu_{\rm exit})\,
                   \rho_{\max}^{\,2}(r_i,\Omega)\,e^{-\tau}\,
                   \mathrm d\Omega,
```
(label: `peierls-class-b-Jn-canonical`)

where the surface-to-observer Jacobian $\mathrm d A_s\,|\mu_s|\,\mathrm d\Omega_{\rm out} = d^{2}\,\mathrm d\Omega_{r_i}$ (with $d = \rho_{\max}$) makes $|\mu_{\rm out}|$ cancel against $|\mu_s|$ so that *every* mode $n$ (including $n = 0$) carries the $\rho_{\max}^{2}$ weight in the observer-centred integrand. The legacy mode-0 `compute_P_esc` *omits* this Jacobian: it returns the unweighted half-sphere outgoing-hemisphere integral $\int_{2\pi^+} e^{-\tau}\,\mathrm d\Omega$ divided by an isotropic-source escape-probability normalisation. The two integrals span **different sub-spaces** of the half-range partial-current basis — $\mu$-weighted (Marshak inner product $\langle f, g\rangle_M = \int_0^1 f g\,\mu\,\mathrm d\mu$) vs unweighted (Lambert inner product $\langle f, g\rangle_L = \int_0^1 f g\,\mathrm d\mu$). The mismatch is exactly the Lambert / Marshak basis change documented for the Class A F.4 closure at `peierls-f4-rank-1-gauge-why` — the algebraic bridge between the two is a non-trivial upper-bidiagonal change-of-basis matrix $M^{(N)}$ (eq. `peierls-M-rank-2`) which becomes a genuine *basis rotation* (not a scalar gauge) at $N \ge 2$. F.4 gets away with the Lambert/Marshak hybrid at rank-1 because the mismatch is the scalar $M^{(1)} = \sqrt{2}/2$ and factors out; the Class B legacy mode-0 / canonical mode-$n \ge 1$ hybrid does the same trick at rank-1, but at rank $\ge 2$ the basis rotation is genuine and the closure is structurally inconsistent.

### Why pure-canonical is not the fix

Switching mode-0 to `compute_P_esc_mode(n=0)` makes the rank-N expansion internally consistent, but the resulting closure is not the right closure: it converges to a wrong limit (~−25 % across all 1R and 2R configurations at high $N$). The reason is that the canonical Marshak partial-current basis with mode-0 weighted by $\rho_{\max}^{2}$ does not reduce to the production rank-1 Mark closure at $N = 1$ — it gives a different (worse) single-mode closure. Pure-canonical breaks the rank-1 regression gate and converges to the wrong answer. **The production rank-1 Mark closure and the canonical rank-N Marshak closure disagree at $N = 1$**; reconciling them requires either re-deriving rank-1 Mark in the canonical basis (which is not the shipped Mark) or re-deriving mode-$n \ge 1$ so that the $N = 1$ truncation is the shipped Mark, with mode-$n \ge 1$ corrections living in a basis where they compose consistently with that mode-0.

</details>

---

## Follow-up directions (2026-04-25 parallel-dispatch investigation)

Three candidate paths to address the Mark uniformity overshoot were investigated by parallel numerics-investigator dispatch.

<details><summary>1. Davison method-of-images (ABANDONED, structurally impossible)</summary>

For sphere with vacuum at center (Davison $u(0)=0$) and white BC at outer surface, the image series

```math
K_{\rm white}(r, r') = \sum_{n} (-1)^{|n|}
    \bigl[E_1(\tau|r - 2nR - r'|) - E_1(\tau|r - 2nR + r'|)\bigr]
```

converges fast (saturates at ~5 image terms, exponentially decaying) but to the **wrong eigenvalue** (−53 % off cp_sphere $k_\infty$ for 1G/1R fuel A). Root cause: method-of-images requires the BC to act **pointwise on the angular flux**. Vacuum ($\psi^- = 0$) and specular ($\psi^-(\Omega) = \psi^+(\Omega - 2(\Omega \cdot \mathbf{n})\mathbf{n})$) qualify. White BC (Mark) does NOT — it re-emits returning current with the AVERAGED angular distribution $\psi^-(\Omega) = J^+/\pi$, which cannot be reproduced by mirror sources. The image series solves the SPECULAR-reflection problem instead. **No method-of-images formulation exists for white BC**, even on a homogeneous sphere; this is a hard structural barrier.

Probes: `derivations/diagnostics/diag_sphere_davison_image_{01..04}_*.py`; agent memory: `.claude/agent-memory/numerics-investigator/issue_132_davison_image_series.md`.

</details>

<details><summary>2. Augmented Nyström (ABANDONED, algebraically equivalent)</summary>

Adding the surface partial current $J^+(\mu)$ as $M$ extra unknowns in an $(N+M)\times(N+M)$ block system, enforcing $J^- = J^+$ as constraint equations rather than a closure approximation. **Verified algebraically equivalent** to the existing $K_{\rm bc} = G \cdot R \cdot P$ Schur reduction in `BoundaryClosureOperator` to machine precision ($6.66\times10^{-16}$ at $M = 1$). The block formulation pre-elimination IS the Schur reduction; eliminating $J^+$ from the bottom block recovers the existing closure with $R = (I - W_{\rm eff})^{-1}$. **The augmented direction provides zero new physics beyond the existing rank-M closure** — the only knob is the choice of $M$-mode basis, and:

- µ-Nyström-collocation is structurally non-convergent ($r$-dependent endpoint singularity at $\mu_{\min}(r) = \sqrt{1 - (r/R)^2}$).
- Marshak Legendre rank-$M$ plateaus at +2.4 % — the same Issue #100 mode-0/mode-$n \ge 1$ normalisation floor that prevented the rank-N Marshak path from being shipped originally.

Probes: `derivations/diagnostics/diag_sphere_augmented_nystrom_{a..e}_*.py`; memory: `.claude/agent-memory/numerics-investigator/issue_132_augmented_nystrom.md`.

</details>

<details><summary>3. Cylinder Hébert + Issue #112 Phase C (RESOLVED for homogeneous)</summary>

The cylinder analog of `compute_P_ss_sphere` derives cleanly:

```math
P_{ss}^{\rm cyl}(\Sigma_t, R) = \frac{4}{\pi}\!\int_0^{\pi/2}
    \cos\alpha\;\mathrm{Ki}_3\!\bigl(2\Sigma_t R\cos\alpha\bigr)\,d\alpha
```
(label: `peierls-cyl-Pss-derivation`)

with $\mathrm{Ki}_3$ arising from analytical integration over the polar angle $\beta$ from the cylinder axis. Multi-region extension mirrors the sphere chord-piecewise integration with the 2-D $h = R\sin\alpha$ impact-parameter geometry. Verified to $<5\times10^{-3}$ against independent Monte Carlo. Shipped as `compute_P_ss_cylinder` with 16 foundation tests at `tests/cp/test_cylinder_pss.py`.

#### Issue #112 Phase C — corrected 3-D G_bc for cylinder

The cylinder rank-1 Mark closure historically used a surface-centric `Ki_1(τ)/d` form lacking the Lambertian projection factor $(R - r\cos\phi)/d$. Row-sum probe quantified the bias:

| Geometry | $\min(K\cdot 1/\Sigma_t)$ | mean | max |
|---|---|---|---|
| Sphere 1G/1R (Hébert) | 0.9993 | 0.9993 | 0.9994 |
| Cylinder 1G/1R (Hébert + buggy G_bc) | 0.8886 | **0.8924** | 0.9089 |
| Cylinder 1G/1R (Hébert + corrected G_bc) | 0.9994 | **0.9996** | 0.9997 |

The corrected 3-D form, derived via SymPy in `derivations/peierls_cylinder_g_bc_3d_derivation.py`, is the observer-centric integral

```math
G_{\rm bc}^{\rm cyl}(r) = \frac{4}{\pi}\!\int_0^\pi
    \mathrm{Ki}_2\!\bigl(\Sigma_t\,d_{\rm 2D}(r, \psi)\bigr)\,d\psi
```
(label: `peierls-cyl-Gbc-3d-final`)

with $d_{\rm 2D}(r, \psi) = -r\cos\psi + \sqrt{R^2 - r^2\sin^2\psi}$ the in-plane backward chord. The $\mathrm{Ki}_2$ arises from analytical integration over the polar angle from the cylinder axis (Knyazev $\mathrm{Ki}_{2+k}$ expansion at $k = 0$).

Shipped as `compute_G_bc_cylinder_3d` and wired into `boundary="white_hebert"` for cylinder (the `NotImplementedError` is now lifted). Cylinder Class B convergence results at BASE quadrature:

| Configuration | cp_cylinder $k_\infty$ | Hébert $k_{\rm eff}$ | err |
|---|---|---|---|
| cyl 1G/1R | 1.500 | 1.4985 | **−0.097 %** |
| cyl 1G/2R *(fuel-A inner / mod-B outer)* | 0.990 | 1.0997 | +11.08 % |
| cyl 2G/1R | 1.875 | 1.8651 | **−0.529 %** |
| cyl 2G/2R | 0.740 | 1.1067 | +49.6 % |

**Same convergence pattern as sphere**: homogeneous configs (1G/1R, 2G/1R) reach <1 % L1 tolerance; heterogeneous configs (1G/2R, 2G/2R) retain the Mark uniformity overshoot. Cylinder 2G/2R is more sensitive than sphere 2G/2R because the cylinder eigenvector is more localized than the sphere eigenvector for the same fuel-mod arrangement. Resolution remains the open question (Sanchez 1977 NSE 64 — see synthesis below).

Probes: `derivations/diagnostics/diag_cylinder_hebert_{pss,keff,diagnose_residual}.py`, `diag_cylinder_g_bc_3d_patched_test.py`; derivation: `derivations/peierls_cylinder_g_bc_3d_derivation.py`; memory: `.claude/agent-memory/numerics-investigator/issue_132_cylinder_hebert.md`.

</details>

### Synthesis

For Class B sphere AND cylinder (after Issue #112 Phase C), the chi-dependent Mark uniformity overshoot is the **intrinsic structural limit** of the rank-N interface-current closure family. Three independent investigations confirm this:

- **Davison method-of-images** is structurally impossible (white BC acts on the *averaged* angular distribution, not pointwise).
- **Augmented Nyström** is algebraically equivalent to the existing rank-M Schur reduction (verified to machine precision at $M=1$).
- **Sanchez 1977 NSE 64 rank-N IC formulation** (PDF read 2026-04-25 — see `.claude/agent-memory/literature-researcher/sanchez_1977_nse64_canonical_ic.md`) uses 3 modes per face = same algebraic class as rank-M Schur with $M=3$ per face. The paper is fixed-source only — its $(I - A\cdot P_{ss})^{-1}$ is *multi-cell coupling*, not a multi-collision $K_\infty$ closure as the earlier ORPHEUS narration asserted. The Hébert $(1 - P_{ss})^{-1}$ is a *different object* (single-cell self-multiplication), not the rank-0 collapse of Sanchez 1977. The empirical bound from Sanchez's Tables VI–VII on heterogeneous LWR cells is ~1–3 % residual with rank-3 — an order of magnitude smaller than the +50 % cyl 2G/2R overshoot, so rank-3 cannot close the gap.

### Live candidate paths (none yet PDF-verified to actually solve the heterogeneous-Mark problem)

To address the overshoot at L1 tolerance for source-localised spectra requires a fundamentally different approach.

- **Sanchez 2002 NSE 142** — double-PN extension of the rank-N IC formulation. Genuinely distinct surface response if the double-PN basis spans more than the per-face piecewise-uniform basis of Sanchez 1977.
- **Bogado Leite 1998 ANE** — orphaned reference per `.claude/agent-memory/literature-researcher/rank_n_ic_curvilinear_literature_leads.md`, candidate for a non-rank-N approach.
- Direct Sn or MOC angular-flux discretisation of the surface reflection — departs from the Peierls integral-equation paradigm.
- Acceptance of the chi-dependent overshoot as the Mark-closure documented limit; restrict L1-tolerance verification to configurations that route through near-uniform-$\sigma_t$ groups (default 2G XS chi=[1, 0]).

To fully resolve the 1G/2R case requires either:

(a) An angular-distribution-preserving closure — the rank-N Marshak path was falsified (Issue #132, this session); rank-N does NOT converge structurally.
(b) The Davison sphere kernel via method-of-images, absorbing the surface reflection analytically into the volume kernel itself via the spherical inversion $r' \mapsto R^2/r'$. This requires literature confirmation that a closed form exists; the method is feasible in principle but no canonical textbook reference was found in the 2026-04-25 literature pull.
(c) Augmented Nyström — adding the surface partial current $J^+$ as an extra unknown in the matrix system, enforcing $J^- = J^+$ as a constraint equation rather than a closure approximation. Issue #100 original suggestion; engineering work but no novel mathematics.

Issue #132 stays OPEN for the 1G/2R limitation; the Hébert closure ships as a substantial partial fix.

---

## Hand-off for the next agent

If you are picking up Class B rank-N MR×MG work:

1. **Read this comment in full** — it is the load-bearing record of why the obvious paths (Davison, Augmented Nyström, pure-canonical Probe G fix) are blocked.
2. **The Hébert (1−P_ss)⁻¹ closure** ships as the substantial partial fix for sphere only. Cylinder + slab raise NotImplementedError. The +10 % sphere 1G/2R overshoot is the documented Mark-closure limitation, pinned by `test_class_b_sphere_hebert_heterogeneous_overshoot_known`. See surviving §`peierls-class-b-sphere-hebert`.
3. **The XFAIL-strict catastrophe pins** at `tests/derivations/test_peierls_rank_n_class_b_mr_mg.py` — `test_class_b_mr_catastrophe_sphere_1g_2r_rank2` (+57 %) and `test_class_b_mr_catastrophe_cylinder_1g_2r_rank2` (+18 %) — flip to **unexpected-pass** when a corrective re-derivation lands. Do not delete those tests; they are the regression gate.
4. **ERR-030** in `tests/l0_error_catalog.md` is the bug catalog entry: failure-mode classification, how it hid (single-region single-group passing rates are degenerate evidence), the lesson (any "rank-N converges" claim must be verified at MR×MG with a non-trivial $\Sigma_t$ breakpoint).
5. **Lessons that generalise (recorded in ERR-030):**
   - The conservation row-sum identity is **not a sufficient gate** when $\Sigma_t$ is uniform (the identity collapses to a tautology); it must be tested with piecewise $\Sigma_t$ to discriminate real conservation from algebraic self-consistency.
   - L19 (signed-error stability under quadrature refinement) does not cover MR/MG configuration refinement — *configuration refinement* is an additional necessary axis on Class B.

References:
- Sanchez & McCormick (1982), *Nuclear Science and Engineering* 80(4), 481–535. §III.F.1, Eqs. 165–169.
- Knyazev (1993), *Atomic Energy* 74(5), DOI 10.1007/BF00844623.
- Sanchez (1977), *Nuclear Science and Engineering* 64.
- Sanchez (2002), *Nuclear Science and Engineering* 142.
- Hébert (2009/2020), *Applied Reactor Physics* (3rd ed.), Ch. 3 §3.8.5, Eq. 3.323.
- Bogado Leite (1998), *Annals of Nuclear Energy*.
- Davison (1957), *Neutron Transport Theory*.

Closing posture: Issue #132 stays OPEN. Hébert (1−P_ss)⁻¹ ships as the sphere partial-fix; cylinder Hébert primitives are derived and shipped (`compute_P_ss_cylinder`, `compute_G_bc_cylinder_3d`) but heterogeneous configs retain the Mark uniformity overshoot. Resolution awaits either Sanchez 2002 double-PN, Bogado Leite 1998, or a fundamentally different closure family.
