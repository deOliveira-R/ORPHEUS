<!-- Posted as part of the post-#138 Peierls documentation cleanup, 2026-04-30 -->
<!-- Source: docs/theory/peierls_nystrom.rst lines 1355-2289 (~935 LoC) -->
<!-- Cleanup commit: [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits) -->

## Issue #117 close-out — Moment-form Nyström archived

**Summary (2026-04-30):**

1. The **moment-form Nyström architecture** (slab + curvilinear) was archived in favor of the unified adaptive `mpmath.quad` path because the verification side of the CP module no longer needs a fast slab K assembly — the single primitive `K_vol_element_adaptive` covers all geometries, and the moment form's performance advantage is only material for a future production discrete-CP solver.
2. The **τ-coordinate Gauss–Laguerre polar form** (the moment-form's own predecessor) was independently archived earlier because Gauss–Laguerre converges only algebraically (`O(n⁻²)`) on integrands polynomial in `e⁻ᵛ` rather than polynomial in `v`, plateauing at relative error `~9×10⁻⁴` regardless of node count and exhibiting non-monotonic GL convergence (smoking-gun of an outer-integrand kink from the dyadic-τ-cap inner subdivision).
3. Current production: slab `K` ships via `_pg.solve_peierls_*` through the unified adaptive `mpmath.quad` path; the moment-form code lives at `derivations/archive/peierls_moments.py` and `peierls_slab_moments_assembly.py` for future reuse; the τ-Laguerre diagnostic scripts remain at `derivations/diagnostics/diag_slab_polar_outer_mu_structure.py` and `diag_slab_polar_glaguerre.py`.

**Surviving Sphinx anchors:** `:label: theory-peierls-moment-form` and `:label: theory-peierls-moment-form-failed-polar` survive as a stub in `docs/theory/peierls_nystrom.rst` pointing here. The math, derivations, conditioning analysis, and numerical evidence below are the architectural reference for any future higher-order discrete-CP solver.

---

<details><summary>Full investigation record (relocated from docs/theory/peierls_nystrom.rst:1355–2289)</summary>

### The unified moment-form statement

For an observer at radial coordinate $r_i$, the unified-form K matrix entry is

```math
K[i, j] = \Sigma_t(r_i)\,C_d
\int_{\Omega_d} W_\Omega(\Omega)\,d\Omega
\int_0^{\rho_{\max}(r_i,\Omega)}
    \frac{1}{\Sigma_t(r'(s,\Omega,r_i))}\,
    \kappa_d\!\bigl(\tau(s)\bigr)\,L_j\!\bigl(r'(s,\Omega,r_i)\bigr)\,ds
```

with $C_d$ the geometry prefactor, $\kappa_d \in \{E_1, \mathrm{Ki}_1, e^{-u}\}$, and $\tau(s) = \int_0^{s}\Sigma_t(r'(t,\Omega,r_i))\,dt$ the cumulative optical depth from the observer.

Substitute the **optical-depth coordinate** $u = \tau(s)$. On each homogeneous segment of the ray $u$ is linear in $s$ with slope $\Sigma_t$, so $ds = du/\Sigma_t$. The inner integral over a single segment with optical-depth range $[u_a, u_b]$ becomes

```math
\int_{u_a}^{u_b}\,\kappa_d(u)\,L_j\!\bigl(r'(u)\bigr)\,
\frac{du}{\Sigma_{t,\text{seg}}^{\,2}},
```

after one factor of $1/\Sigma_t$ from the $1/\Sigma_t(r'_{ikm})$ already in the K-source form and one from $ds = du/\Sigma_{t,\text{seg}}$. (The unified-form prefactor $\Sigma_t(r_i)$ absorbs the second $1/\Sigma_t$ so the matrix entry stays dimensionless.)

If the panel basis $L_j$ is polynomial of degree $p-1$ in $u$, expand it in monomials,

```math
L_j(r'(u)) = \sum_{m=0}^{p-1} c^{(j)}_m\,u^m,
```

so the segment integral collapses to the inner product

```math
\boxed{\;
\int_{u_a}^{u_b}\!\kappa_d(u)\,L_j(r'(u))\,du
=
\sum_{m=0}^{p-1} c^{(j)}_m\,\bigl[J_m^{\kappa_d}(u_b) - J_m^{\kappa_d}(u_a)\bigr]
=
\langle\,\mathbf c^{(j)}\,,\,\mathbf M^{(\kappa_d)}\,\rangle.
\;}
```

The two factors play distinct roles:

- $\mathbf c^{(j)}$ — **basis coefficients**, geometry-specific through the chord map $r'(u)$ but kernel-independent. Built once per (observer, segment) pair.
- $\mathbf M^{(\kappa_d)}_m = J_m^{\kappa_d}(u_b) - J_m^{\kappa_d}(u_a)$ — **kernel moments**, geometry-independent functions of the optical-depth endpoints $(u_a, u_b)$ only. The cumulative moment $J_m^{\kappa_d}(z)\equiv\int_0^z u^m\,\kappa_d(u)\,du$ admits a closed-form recursion for every kernel of interest.

K-matrix assembly pseudocode:

```text
for each observer i:
    for each source panel p:
        build segments of the ray from x_i intersecting panel p
        for each segment with [u_a, u_b]:
            M       = kernel_moments(u_a, u_b, p_order)         # closed form
            c[a,m]  = basis_coeffs_in_u(panel_nodes, segment)   # Vandermonde solve
            K[i, panel_node_a] += prefactor · c[a,:] · M
```

with no inner quadrature.

When the chord $r'(u)$ is non-polynomial in $u$ (cylinder, sphere), $L_j(r'(u))$ is no longer a polynomial in $u$, and the closed-form contraction requires *interpolation* to recover polynomial coefficients. This recovers the standard Gauss–Legendre Nyström — but applied at the kernel-moment-weighted level rather than the kernel-evaluation level.

### Closed-form moment recursions

The cumulative kernel moments $J_m^{\kappa_d}(z)$ are derived once per kernel by integration by parts. The implementations live in `orpheus.derivations.peierls_moments` (now archived) with element-wise gates against `mpmath.quad` to $10^{-15}$ relative.

**Slab — $E_1$ moments.** From $E_1'(u) = -e^{-u}/u$ and $(u^{m+1})' = (m+1)\,u^m$, integration by parts gives

```math
J_m^{E_1}(z) \equiv \int_0^z u^m\,E_1(u)\,du
= \frac{z^{m+1}\,E_1(z) + \gamma(m+1, z)}{m+1},
```

where $\gamma(a, z) = \int_0^z t^{a-1} e^{-t}\,dt$ is the lower incomplete gamma function. The intermediate step:

```math
\int_0^z u^m E_1(u)\,du
= \left[\frac{u^{m+1}}{m+1}\,E_1(u)\right]_0^z
  - \frac{1}{m+1}\int_0^z u^{m+1}\cdot\!\left(-\frac{e^{-u}}{u}\right)\!du
= \frac{z^{m+1} E_1(z)}{m+1} + \frac{1}{m+1}\int_0^z u^m e^{-u}\,du.
```

The boundary term at $u=0$ vanishes because $u^{m+1}E_1(u) \sim -u^{m+1}\ln u \to 0$ for $m \ge 0$. This identity is LewisMiller1984 Appendix C, Hebert2020 §3.2–3.3, and AbramowitzStegun1964 §5.1.32.

**Cylinder — $\mathrm{Ki}_1$ moments.** From $\mathrm{Ki}_n'(u) = -\mathrm{Ki}_{n-1}(u)$ and $\int\mathrm{Ki}_n(u)\,du = -\mathrm{Ki}_{n+1}(u)$, repeated integration by parts of $\int_0^z u^m \mathrm{Ki}_1(u)\,du$ gives

```math
J_m^{\mathrm{Ki}_1}(z) = -\sum_{q=0}^{m}
  \frac{m!}{(m-q)!}\,z^{m-q}\,\mathrm{Ki}_{q+2}(z)
  + m!\,\mathrm{Ki}_{m+2}(0).
```

Each integration by parts trades one power of $u$ for one increment of the Bickley index, telescoping $u^m\mathrm{Ki}_1$ all the way down to $m!\,\mathrm{Ki}_{m+2}$ in the boundary terms. The endpoint constants $\mathrm{Ki}_n(0) = (\sqrt{\pi}/2)\,\Gamma(n/2)/\Gamma((n+1)/2)$ come from Wallis' formula evaluated at the explicit definition $\mathrm{Ki}_n(0) = \int_0^{\pi/2} \cos^{n-1}\theta\,d\theta$. The recursion is implicit in Stamm1983 Chapters 4–6 and is restated in modern notation by Hebert2020 §3.4–3.5; underlying Bickley identities in AbramowitzStegun1964 §11.2 and Bickley.

The $\mathrm{Ki}_n(0)$ constants are computed at the caller's mpmath `workdps` precision via the Wallis closed form because the small-$z$ regime exhibits cancellation between the boundary $m!\,\mathrm{Ki}_{m+2}(0)$ term and the sum of $u^{m-q}\mathrm{Ki}_{q+2}(z)$ terms, amplifying any roundoff in the boundary constants. A float-precision constant would limit the moment recursion to $\sim10^{-9}$ relative at small $z$; the mpmath-native constant gives full `dps` precision.

**Sphere — exponential moments.** $\int_0^z u^m e^{-u}\,du = \gamma(m+1, z)$ directly; no integration by parts needed.

Per-segment differences (used by the K assembly):

```math
\mathbf M^{(\kappa)}_m \equiv
\int_{u_a}^{u_b}\!u^m\,\kappa(u)\,du
= J_m^{\kappa}(u_b) - J_m^{\kappa}(u_a),
\qquad m = 0, 1, \dots, p-1.
```

For the slab, computing $\mathbf M^{(E_1)}$ requires exactly two $E_1$ evaluations and $2p$ lower-incomplete-gamma evaluations — entirely closed-form, no quadrature.

### Slab specialisation (production at the time)

The slab is the "easy" geometry because the chord map is linear: $r'(s) = x_i + s\,\sigma$, where $\sigma = \pm 1$ is the ray-direction sign. With $u = \Sigma_{t,\text{seg}}\,s$, $r'(u) = x_i + \sigma\,u/\Sigma_{t,\text{seg}}$ is **linear in** $u$. A panel Lagrange basis $L_a(r')$ of degree $p-1$ in $r'$ remains degree $p-1$ in $u$.

Implementation `_build_volume_kernel_slab_moments` (now archived) processes each (observer, source panel) pair as follows:

**Step 1 — Walk the ray.** For each observer $x_i$, walk the ray to the source panel and accumulate the cumulative-$\tau$ offset $\Delta = \Sigma_{t}\,(x_l - x_i)$ (right-going) or $\Sigma_{t}\,(x_i - x_r)$ (left-going), using a piecewise-linear cumulative-tau table over the material regions of the slab:

```math
\Delta = \int_{x_i}^{\text{ray-entry into panel}} \Sigma_t(\xi)\,d\xi.
```

This handles heterogeneous material boundaries between observer and source panel without unnecessary subdivision. (Source panels do not span material boundaries by construction.)

**Step 2 — Self-panel splitting.** When the observer sits inside the source panel ($x_l \le x_i \le x_r$), the ray-from-observer leaves the panel in two opposite directions. Split into two ray pieces:

```math
\text{piece A: } [x_i, x_r],\ \sigma=+1,\ \Delta = 0; \qquad
\text{piece B: } [x_l, x_i],\ \sigma=-1,\ \Delta = 0.
```

Each piece is $C^\infty$ in $u$; no log-singularity subtraction is needed because $E_1(u)$ itself is $\sim -\ln u + O(1)$ near the origin and our **moment** recursion already absorbs that singularity into the closed form (the boundary term $z^{m+1}E_1(z)/(m+1)$ carries the log analytically). This is the central reason the slab-moment form succeeds where the legacy $E_1$ Nyström needed Atkinson's product-integration recipe (Atkinson1997): the closed form integrates the singularity in *symbolic* form, so the residual integrand the quadrature would have seen is already exactly $\gamma(m+1, z) / (m+1)$ — finite and smooth.

**Step 3 — Build the per-segment moment vector.** With segment optical-depth range $[u_a = \Delta,\,u_b = \Delta + \Sigma_{t,\text{panel}}\,(x_r - x_l)]$ (or piece-specific endpoints for self-panel halves), compute

```math
\mathbf M_m = J_m^{E_1}(u_b) - J_m^{E_1}(u_a),
\qquad m = 0, 1, \dots, p-1,
```

via the closed form evaluated at mpmath `dps` precision (default 30).

**Step 4 — Vandermonde solve for cardinal Nyström weights.** The panel nodes $\xi_a$ define cardinal Lagrange polynomials $L_a(r') = \prod_{b\neq a}(r' - \xi_b)/(\xi_a - \xi_b)$. In the segment's optical-depth coordinate, each node sits at

```math
u_k = \Delta + \Sigma_{t,\text{panel}}\,\sigma\,(\xi_k - x_{\text{a,seg}})
```

(with $x_{\text{a,seg}} = x_l$ for $\sigma=+1$, $x_{\text{a,seg}} = x_r$ for $\sigma=-1$). Solve

```math
\mathbf V^\top\,\mathbf w = \mathbf M,
\qquad \mathbf V[k, m] = u_k^m,
```

where $\mathbf V$ is the $p\times p$ Vandermonde matrix at the panel-node $u$-coordinates. The transpose comes from the dual cardinal-vs-monomial relation: $L_a(r'(u)) = \sum_m c^{(a)}_m\,u^m$ has coefficients $c^{(a)} = (\mathbf V^{-1})_{a,:}$, so $w_a = \sum_m (\mathbf V^{-1})_{a,m}\,M_m = (\mathbf V^{-\top}\mathbf M)_a$.

The system is solved at mpmath `dps` precision with `mpmath.lu_solve`. The Vandermonde matrix is notoriously ill-conditioned in floating point — for distant source panels (large $\Delta$) the panel-node $u$-coordinates cluster (relative spread $(\xi_p - \xi_0)\,\Sigma_t / \Delta$ becomes small), and a float-precision solve loses 6-10 digits. Solving at $\mathrm{dps}=30$ gives single-digit ULP loss in the final K entry — well below the $10^{-12}$ test gate.

**Step 5 — Accumulate into K.** With Nyström weights in hand:

```math
K[i, j_{\text{start}+a}] \mathrel{+}=
\frac{\Sigma_t(r_i)\,C_d}{\Sigma_{t,\text{panel}}}\;w_a,
\qquad a = 0, 1, \dots, p-1,
```

where $C_d = 1/2$ for the slab. The full algorithm is $O(N \cdot N_{\text{panels}} \cdot p^3)$ — linear in (observer, source-panel) pairs and cubic in panel order through the LU. For typical $N=24, p=4$ it builds in $\sim 60$ ms at $\mathrm{dps}=30$; at $p=6, N=24$, $\sim 200$ ms.

The slab moment-form is **exact in the inner integration**: no inner quadrature, no inner discretization error, no inner convergence sweep. Sources of error: (a) panel basis order $p-1$ (truncation, converges spectrally for smooth solutions), (b) mpmath precision of moment evaluations and Vandermonde solve (default $\mathrm{dps}=30$ gives $\sim10^{-28}$ relative).

</details>

<details><summary>Falsified-experiment-of-falsified-experiment: τ-Laguerre polar-form (the moment form's own predecessor)</summary>

### Why the predecessor τ-Laguerre polar form failed

The first attempt at unifying the slab into the polar-form Nyström framework (commit 4395cb8, "feat(derivations): add slab-polar as first-class CurvilinearGeometry kind") used a **τ-coordinate Gauss–Laguerre** outer integration in the substitution $v = -\ln|\mu|$. **That implementation was retired** because it does not converge to machine precision — it plateaus at relative error $\sim 9 \times 10^{-4}$ regardless of how many quadrature nodes are added. Diagnostic scripts `derivations/diagnostics/diag_slab_polar_outer_mu_structure.py` and `diag_slab_polar_glaguerre.py` isolate the cause.

**The slab-polar formulation.** With observer at $x_i$ and direction-cosine $\mu$, the polar-form slab K is

```math
K[i, j] = \tfrac{1}{2}\,\Sigma_t(r_i)
\int_{-1}^{1}\!d\mu
  \int_0^{\rho_{\max}(x_i,\mu)}
    \frac{e^{-\Sigma_t \rho}}{\Sigma_t(r')}\,L_j(x_i + \rho\mu)\,d\rho.
```

Substitute $v = -\ln|\mu|$ so $\mu = \pm e^{-v}$ with Jacobian $|d\mu/dv| = e^{-v}$. The outer integral on each branch becomes $\int_0^\infty g(e^{-v})\,e^{-v}\,dv$, where $g(\mu)$ is the inner ρ-integral as a function of $\mu$. Apply Gauss–Laguerre to the outer $v$.

**The defect.** Gauss–Laguerre is *spectrally* accurate when the integrand on $(0, \infty)$ has the form $p(v)\cdot e^{-v}$ with $p$ polynomial in $v$ — that is the rule's definition of "polynomial-exact." Our integrand has the form

```math
g(e^{-v})\,e^{-v},
```

where $g(\mu)$ admits the small-$\mu$ Laplace expansion

```math
g(\mu) = \frac{L_j(x_i)}{\Sigma_t}
       + \mu\,\frac{L_j'(x_i)}{\Sigma_t^2}
       + \mu^2\,\frac{L_j''(x_i)}{\Sigma_t^3}
       + \dots
```

(from differentiating the absorption integral by parts in $\mu$), so $g(\mu)$ is a polynomial of degree $p-1$ *in* $\mu = e^{-v}$. The full integrand is therefore a sum of terms $e^{-kv}$ for $k = 1, 2, \dots, p$, **not** a polynomial in $v$.

Gauss–Laguerre $n$-point integrates $p(v)\,e^{-v}$ exactly for $\deg p \le 2n - 1$; on a generic $e^{-kv}$ integrand it converges only **algebraically** at rate $O(n^{-2})$ from the remainder bound for non-polynomial integrands (Gauss-Laguerre is asymptotic in $1/n$ for any $e^{-kv}$ with $k\neq 1$). Our integrand mixes $k=1$ through $k=p$, so the leading-order error is set by the worst $k$.

**Diagnostic confirmation** (table from `diag_slab_polar_outer_mu_structure.py`, $L=1$, $\Sigma_t=1$, two-panel $p=4$, exact-inner `mpmath.quad` so the outer is isolated). Relative error in $K[0,0]$ against the adaptive reference $K[0,0] = 1.107\times10^{-1}$:

| $n_v$ | v-Laguerre | v-GL on $[0, 30]$ | GL on $[-1,0]\cup[0,1]$ |
|------:|-----------:|------------------:|------------------------:|
|   8   |     6.6e-3 |            5.2e-3 |                  1.1e-3 |
|  16   |     1.6e-3 |            1.3e-3 |                  7.4e-5 |
|  32   |     4.0e-4 |            3.4e-4 |                  6.2e-6 |
|  64   |     1.0e-4 |            8.7e-5 |                  7.1e-4 |

Two readings of this table are essential:

1. The **v-Laguerre and v-GL** rates are both $O(n^{-2})$ — confirming the polynomial-in-$e^{-v}$ defect is substitution-independent (any rule chosen for the polynomial-in-$v$ weight $e^{-v}$ shares the same exactness class).
2. The **standard GL** convergence is *non-monotonic*: it drops to $6\times 10^{-6}$ at $n=32$ then *jumps back up* to $7\times 10^{-4}$ at $n=64$. This is the smoking gun of an outer-integrand kink that the GL node placement straddles differently at each $n$. The kink is the predicted signature of the dyadic-$\tau$-cap on the inner subdivision (hypothesis H3 in the diagnostic docstring) — a discrete change in the number of inner sub-intervals creates a discontinuity in $g(\mu)$ as $\mu$ moves past the cap-trigger threshold. The moment-form sidesteps this entirely because the inner integral is closed-form; no cap is needed.

**Generalised Laguerre does not help.** The diagnostic `diag_slab_polar_glaguerre.py` sweeps the $\alpha$ parameter of the generalised Laguerre weight $v^\alpha\,e^{-v}$ (intuition: shifting nodes toward small $v$ should help because the integrand is concentrated there). Result (relative error in $K[0,0]$, $n_v=32$, exact inner):

| $n_v$ | $\alpha = -0.5$ | $\alpha = 0.0$ | $\alpha = +0.5$ | $\alpha = +1.0$ | $\alpha = +2.0$ | $\alpha = +5.0$ |
|------:|----------------:|---------------:|----------------:|----------------:|----------------:|----------------:|
|   8   |          6.5e-3 |         6.6e-3 |          7.1e-3 |          8.1e-3 |          1.1e-2 |          2.5e-2 |
|  16   |          1.6e-3 |         1.6e-3 |          1.8e-3 |          2.2e-3 |          3.4e-3 |          9.5e-3 |
|  32   |          4.0e-4 |         4.0e-4 |          4.5e-4 |          5.7e-4 |          1.0e-3 |          3.6e-3 |
|  64   |          1.0e-4 |         1.0e-4 |          1.1e-4 |          1.5e-4 |          2.8e-4 |          1.1e-3 |

Every value of $\alpha$ is in the same algebraic-convergence class as standard Laguerre, with $\alpha > 0$ strictly worse (node clustering toward large $v$ is precisely the wrong place for an integrand concentrated at small $v$) and $\alpha < 0$ giving no measurable improvement. The theoretical reason is that **all Laguerre flavours, regardless of** $\alpha$, **are spectral only for polynomial-in-**$v$ **integrands** — the weight $v^\alpha e^{-v}$ does not change the rule's polynomial-exactness class. Our integrand is polynomial in $e^{-v}$, not in $v$. There is no Laguerre flavour that fixes this.

**The architectural realisation.** The legacy $E_1$ Nyström does not exhibit this convergence pathology because it integrates the angular variable **analytically** — the classical slab equation already has $E_1(\Sigma_t|x-x'|)$ as the kernel, and $E_1$ *is* the analytical result of $\int_0^\infty e^{-u}/u\,du$ along the angular direction. The polar form was undoing that analytical integration and trying to recover it numerically.

The unified moment form re-uses $E_1$ *as the kernel* (the way the classical slab Nyström does) but re-writes the angular-then-radial integration as a single optical-depth-coordinate integral with a polynomial source. The result is the polynomial-coefficient × kernel-moment contraction — exact, closed-form, no quadrature.

</details>

<details><summary>Curvilinear status (cylinder and sphere) — kernel-natural moment form</summary>

For the cylinder and sphere geometries the chord map $r'(\rho) = \sqrt{r_i^2 + \rho^2 - 2 r_i \rho \cos\Omega}$ is non-linear in $\rho$, hence non-linear in the optical-depth coordinate $u = \Sigma_t\,\rho$ even within a homogeneous segment. A panel basis $L_j(r')$ polynomial of degree $p-1$ in $r'$ is **not** polynomial in $u$ once composed with the chord map. The closed-form contraction therefore does not give a closed-form K entry; the polynomial coefficients of $L_j(r'(u))$ are not a finite list.

The unified architecture handles this — implemented in `_build_volume_kernel_curvilinear_moments` (archived) — by **interpolating** $L_j(r'(u))$ at $n_\rho$ Gauss–Legendre nodes inside each segment, recovering polynomial coefficients $c^{(j)}_m$ from a Vandermonde solve at those $u$-nodes, and contracting against the closed-form $\mathrm{Ki}_n / e^{-u}$ moments:

```math
K[i, j] \mathrel{+}=
\frac{\Sigma_t(r_i)\,C_d}{\Sigma_{t,\text{seg}}}\,
\sum_{k=0}^{n_\rho-1} w_k\,L_j\!\bigl(r'(u_k)\bigr),
```

with $w_k = \langle (\mathbf V^{-\top})_{k,:}, \mathbf M^{(\kappa)}\rangle$ the Nyström weight at node $u_k$.

**Same architecture, same moment recursions, same Vandermonde solve**, but the polynomial coefficients are obtained by interpolation rather than from the panel basis directly. We call this the *kernel-natural moment form*.

**Why this is not a regression.** The legacy curvilinear path evaluates the kernel $\kappa_d(u_k)$ at each Gauss–Legendre node and quadrature-sums. The kernel evaluation at isolated nodes is just a numerical approximation of the closed-form moment $\int_{u_a}^{u_b}\!u^m \kappa_d(u)\,du$ accurate to the GL polynomial-exactness class $(2 n_\rho - 1)$ in $u$. The natural-kernel moment form **uses the closed form directly** — one architectural step *more* accurate per node, reducing required $n_\rho$ for a given precision. Otherwise the spectral convergence in $n_\rho$ is identical (set by the smoothness of $L_j(r'(u))$).

Architecture instantiation by geometry:

| Geometry | Kernel $\kappa_d$ | Moment recursion | Polynomial coeffs | Inner quadrature |
|----------|-------------------|------------------|-------------------|------------------|
| Slab | $E_1$ | $J_m^{E_1}(z)$ closed form | panel-cardinal (closed) | **none** |
| Cylinder | $\mathrm{Ki}_1$ | $J_m^{\mathrm{Ki}_1}(z)$ closed form | GL-$u$-cardinal | $n_\rho$ GL on $L_j(r'(u))$ |
| Sphere | $e^{-u}$ | $\gamma(m+1, z)$ | GL-$u$-cardinal | $n_\rho$ GL on $L_j(r'(u))$ |

</details>

<details><summary>Numerical evidence (4 gates, all passed at the time of archival)</summary>

The slab moment-form K matrix was gated against three independent references in `tests.derivations.test_peierls_slab_moments` (L1 equivalence, all marked `@pytest.mark.verifies("peierls-equation")`). All gates passed on commit `investigate/peierls-solver-bugs`.

**Gate 1 — moment vs legacy E₁ Nyström, homogeneous slabs.** Element-wise relative error $< 10^{-12}$ across five parametrizations:

| $L$ | $\Sigma_t$ | $n_{\text{panels}}$ | $p$ | status |
|----:|-----------:|--------------------:|----:|--------|
| 1.0 |        1.0 |                   2 |   4 | pass |
| 1.0 |        1.0 |                   4 |   6 | pass |
| 2.0 |        0.5 |                   3 |   4 | pass |
| 1.0 |        5.0 |                   2 |   4 | pass (optically thick, $\Sigma_t L = 5$) |
| 0.5 |        2.0 |                   4 |   4 | pass (short slab) |

The $p=6$ parametrisation exercises the moment recursion at $m=0,1,\dots,5$ and the $6\times 6$ Vandermonde — the hardest to condition in float arithmetic. The mpmath-$\mathrm{dps}=30$ solve preserves $\sim 28$ digits, comfortably below the $10^{-12}$ test gate.

**Gate 2 — moment vs legacy E₁ Nyström, heterogeneous slabs.** Same gate, three two-region parametrizations (panel grid respects material boundary):

| region thicknesses | $\Sigma_t$ regions | $n_{\text{panels per region}}$ | $p$ | status |
|--------------------|--------------------|-------------------------------:|----:|--------|
| [1.0, 1.0]         | [1.0, 0.5]         |                              2 |   4 | pass |
| [0.5, 1.5]         | [2.0, 0.3]         |                              3 |   4 | pass |
| [1.0, 0.5]         | [0.8, 4.0]         |                              4 |   4 | pass |

The third parametrisation is the most demanding: thin optically-thick region ($\Sigma_t L = 4 \cdot 0.5 = 2$) directly adjacent to thicker optically-thin region ($\Sigma_t L = 0.8$). The cumulative-tau walker must compose the two regions correctly; a sign-flip or material-region indexing bug would surface as a per-cross-region $K[i,j]$ discrepancy. All $(i, j)$ entries pass $10^{-12}$.

**Gate 3 — moment vs adaptive polar reference, element-wise.** Spot-check five entries against `slab_polar_K_vol_element` (nested adaptive `mpmath.quad`), $L=1$, $\Sigma_t=1$, 2 panels $p=4$. Tolerance $10^{-10}$:

| $(i, j)$    | role                                          | status |
|-------------|-----------------------------------------------|--------|
| $(0, 0)$    | leftmost-observer self-contribution            | pass |
| $(0, N-1)$  | leftmost observer to rightmost source          | pass |
| $(N/2, N/2)$| middle-of-slab self-contribution               | pass |
| $(N-1, 0)$  | rightmost observer to leftmost source          | pass |
| $(1, 3)$    | cross-panel arbitrary entry                    | pass |

The $10^{-10}$ floor (vs $10^{-12}$ in Gates 1–2) reflects the adaptive `mpmath.quad` reference's own discretisation limit; the moment form itself is closer to mpmath ULP.

**Gate 0 — moment recursion vs mpmath.quad (term verification).** L0 gates in `tests.derivations.test_peierls_moments`, 32 parametrisations ($z \in \{10^{-3}, 10^{-2}, 0.1, 0.5, 1, 2.5, 5, 10, 25\}$, $m \in \{0, 1, \dots, 6\}$):

| moment family | tolerance | basis |
|---------------|----------:|-------|
| $J_m^{E_1}$ (slab) | $10^{-13}$ | integration by parts → closed form |
| $J_m^{\mathrm{Ki}_1}$ (cylinder) | $10^{-12}$ | repeated IBP → closed form |
| $J_m^{e^{-u}}$ (sphere) | $10^{-15}$ | direct via $\gamma(m+1, z)$ |

The slightly looser $\mathrm{Ki}_1$ tolerance reflects the small-$z$ cancellation between the boundary $m!\,\mathrm{Ki}_{m+2}(0)$ term and the sum of $\mathrm{Ki}_{q+2}(z)$ terms; even with mpmath-native $\mathrm{Ki}_n(0)$ constants, the cancellation reduces working precision by 2-3 digits at $z = 10^{-3}$. This does not affect the K-matrix gate because segment endpoint optical depths in any practical slab/cylinder geometry are $> 10^{-3}$ by orders of magnitude.

Together these four gates established that the moment-form K matrix **equals** the legacy adaptive references at the $10^{-12}$ level, with the underlying recursions verified at $10^{-13}$ or better — well below the power-iteration eigenvalue solver tolerance ($10^{-10}$ typical), so the moment form was a drop-in replacement at archive time.

</details>

<details><summary>Performance characteristics</summary>

The moment-form slab K assembly is dominated by two costs:

1. The closed-form moment vector evaluation `slab_segment_moments_mp`: $2(p+1)$ `mpmath.gammainc` evaluations per segment plus two `mpmath.expint(1, ·)` evaluations.
2. The $p \times p$ Vandermonde LU solve at mpmath `dps` precision.

For a representative $N = 24, p = 4, n_{\text{panels}} = 6$ homogeneous problem at $\mathrm{dps} = 30$, K assembly takes $\sim 60$ ms wall-clock; for $p = 6$ the same problem takes $\sim 200$ ms. By comparison the legacy $E_1$ Nyström at the same precision and same K size takes $\sim 800$ ms, dominated by adaptive `mpmath.quad` calls (each adaptive call spends $\sim 50$ integrand evaluations at $\mathrm{dps}=25$). The moment form is therefore both **simpler** (closed form vs adaptive) and **faster** (no adaptive sub-grid management).

A future optimisation: perform the Vandermonde solve in **float precision** when the panel-node $u$-spread is wide enough that the system is well-conditioned (i.e. for self-panels and near-observer source panels). The conditioning is $\kappa_2(\mathbf V) \sim ((\xi_p - \xi_0)\,\Sigma_t / \langle u \rangle)^{-(p-1)}$ roughly; for $\langle u \rangle \lesssim 5$ (within a few mean free paths) the condition number is $\le 10^{8}$ at $p=4$ and a float solve is safe. For more distant source panels the mpmath solve is required; a hybrid switch could shave another $\sim 3\times$ off the wall-clock without compromising precision. Not currently scheduled.

</details>

### References (preserved verbatim)

- LewisMiller1984 Appendix C — slab CP polynomial-source integration-by-parts identity for $\int_0^z u^m E_1(u)\,du$. Closest single-source statement of the slab moment recursion.
- Hebert2020 §3.2-3.5 — modern restatement of slab (§3.2-3.3) and cylindrical (§3.4-3.5) polynomial-source CP recursions in the language of integration-by-parts moment expansions. Hébert is the most accessible textbook reference.
- Stamm1983 Chapters 4-6 — canonical derivation of the cylinder Bickley-recursion family, including the Wallis closed form $\mathrm{Ki}_n(0) = (\sqrt\pi/2)\,\Gamma(n/2)/\Gamma((n+1)/2)$.
- AbramowitzStegun1964 §5.1.32 — underlying recursion $E_n'(u) = -E_{n-1}(u)$ and boundary behaviour $u^{m+1} E_1(u) \to 0$. §11.2 covers the Bickley-Naylor $\mathrm{Ki}_n$ identities.
- Bickley — original 1935 introduction of the Bickley-Naylor family; identities $\mathrm{Ki}_n'=-\mathrm{Ki}_{n-1}$ and $\int\mathrm{Ki}_n\,du = -\mathrm{Ki}_{n+1}$.

### Production decision (2026-04-30)

- **Archived:** `derivations/archive/peierls_moments.py` and `derivations/archive/peierls_slab_moments_assembly.py` retain the implementation for future production discrete-CP work.
- **Current production:** Slab K ships via `_pg.solve_peierls_*` (the unified adaptive `mpmath.quad` path through `K_vol_element_adaptive`). The verification side does not need a fast slab K — adaptive is fast enough at the verification regime ($N < 50$, dps=25).
- **Diagnostic scripts retained:** `derivations/diagnostics/diag_slab_polar_outer_mu_structure.py` and `diag_slab_polar_glaguerre.py` reproduce the τ-Laguerre tables above for any future session that questions the GL convergence claim.
- **Sphinx stub:** `:label: theory-peierls-moment-form` and `:label: theory-peierls-moment-form-failed-polar` survive in `docs/theory/peierls_nystrom.rst` as a 25-LoC pointer to this comment, the archive, and the diagnostic scripts.

### Cross-links

- Surviving Sphinx anchor: `:ref:\`theory-peierls-moment-form\`` (in `docs/theory/peierls_nystrom.rst`)
- Failed-polar anchor: `:ref:\`theory-peierls-moment-form-failed-polar\``
- Code archive: `derivations/archive/peierls_moments.py`, `derivations/archive/peierls_slab_moments_assembly.py`
- Diagnostics: `derivations/diagnostics/diag_slab_polar_outer_mu_structure.py`, `derivations/diagnostics/diag_slab_polar_glaguerre.py`
- Cleanup commit: [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
- Moment-form landing commit: `investigate/peierls-solver-bugs` (2026-04-19)
- τ-Laguerre landing commit: `4395cb8`
