# Quadrature landscape survey — ORPHEUS

**Branch** `refactor/operator-strategy-layers` @ `1659d756`. **Date** 2026-08-01.
**Mode** read-only inventory. No production code or test touched.

`[M]` = I measured it (probe script + output quoted). Everything else is
read-from-source; inference is labelled *inference*. Docstring claims are quoted
as CLAIMS and separately adjudicated (L33: a docstring is not evidence).

---

## 0. Tree state (L-012 discipline)

`git status --short` at open: 11 modified (all `.claude/` memory + `vv-principles`
skill + `tests/sn/_test_helpers.py`), plus untracked scratch/ and:

- `orpheus/numerics/roots_of_unity.py` — **UNTRACKED, in-flight** (the #325
  algebraic-node work). Treated below as a *proposed* realization tagged
  `(untracked, in-flight)`, NOT as landed code.
- `tests/numerics/test_roots_of_unity.py`,
  `tests/sn/sweep/curvilinear/test_alpha_closed_form.py`,
  `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py`,
  `tests/sn/verification/mms/test_mms_ordering_blindness.py` — untracked, in-flight
  (#325/#326).

Everything else cited is at HEAD.

---

## 1. HEADLINE — the three findings that change the decision

### F1 `[M]` `level_symmetric_sn` **over-claims its degree of exactness by up to 12**, and the registry selector consumes the false number

The tag is `degree_of_exactness = sn_order - 1`
(`orpheus/numerics/quadrature/rules_sphere.py:278`). Measured actual degree
(largest `d` such that EVERY monomial `x^a y^b z^c`, `a+b+c ≤ d`, integrates to
`< 1e-11` relative error):

```
  N=  2  n_nodes=   8  claimed=  1  MEASURED=  3   UNDER-CLAIM
  N=  4  n_nodes=  24  claimed=  3  MEASURED=  3   OK
  N=  6  n_nodes=  48  claimed=  5  MEASURED=  3   *** OVER-CLAIM ***
  N=  8  n_nodes=  80  claimed=  7  MEASURED=  3   *** OVER-CLAIM ***
  N= 10  n_nodes= 120  claimed=  9  MEASURED=  3   *** OVER-CLAIM ***
  N= 12  n_nodes= 168  claimed= 11  MEASURED=  3   *** OVER-CLAIM ***
  N= 16  n_nodes= 288  claimed= 15  MEASURED=  3   *** OVER-CLAIM ***
```

The measured degree is **3 for every N**, and it does not improve with N.
Root cause, read from the implementation: the weights are
**equal within an octant** —

```python
# orpheus/numerics/quadrature/rules_sphere.py:188
w_octant = 4.0 * np.pi / (8.0 * n_octant)
```

`[M] n_distinct_weights == 1` for every N ∈ {2,4,6,8,16}. That is *not* the
Carlson–Lathrop construction the docstring cites; real LQ_N solves moment
conditions for a **set of distinct point weights**. The even-moment error:

```
  N=  2  mu^2 err=0.0e+00   mu^4 err=4.444e-01   mu^6 err=7.407e-01
  N=  4  mu^2 err=2.1e-16   mu^4 err=1.667e-01   mu^6 err=2.870e-01
  N=  6  mu^2 err=0.0e+00   mu^4 err=5.382e-02   mu^6 err=9.816e-02
  N=  8  mu^2 err=0.0e+00   mu^4 err=4.306e-02   mu^6 err=8.327e-02
  N= 16  mu^2 err=0.0e+00   mu^4 err=8.000e-02   mu^6 err=1.514e-01
```

`Σ w μ⁴` is **8 % wrong at S₁₆ and non-monotone in N** (5.4 % at S₆ → 4.3 % at
S₈ → 8.0 % at S₁₆). This rule does not converge in the polynomial-exactness sense
at all.

And the degree-3 that *is* achieved is **free from the symmetry alone**, not from
the level construction `[M]`:

```
  single RANDOM O_h orbit (48 pts, equal weights, sum=4pi):
     degree 0..3 worst rel err <= 1.3e-15;  degree 4 worst rel err = 5.41e-01
```

Any `O_h`-invariant node set with weights summing to `4π` is degree-3 exact
(on `S²`, `x²+y²+z²≡1` plus `O_h` forces `Σwx² = Σwy² = Σwz² = 4π/3`). So the
`level_symmetric` tag currently reports a number that (a) is wrong, and (b) even
when right at N=4 reports only what `invariance_group=OctahedralOh` already
implies.

**The false number has a live consumer** `[M]`. `select_quadrature` inverts the
formula (`_ls_sn_invert`, `registry.py:355`) and does choose this rule:

```
  cartesian2d d= 3 level_structured -> LevelSymmetricSN({'sn_order': 4})
  cartesian2d d= 7 level_structured -> LevelSymmetricSN({'sn_order': 8})
  cartesian2d d=15 level_structured -> LevelSymmetricSN({'sn_order': 16})
```

The last line: the user asked for degree 15, the selector returned a rule whose
measured degree is 3, and the returned measure's `degree_of_exactness` attribute
says 15.

**Why no test catches it (vv Mode 7 — fixture-restricted assertion).**
`tests/numerics/test_rules_sphere.py:111` asserts only that the *tag* equals
`sn_order - 1` — it pins the label, not the property. The only exactness tests in
that file are `test_level_symmetric_integrates_constants` (degree 0) and
`test_level_symmetric_quadratic_isotropy` (degree 2) — both inside the sector the
symmetry makes free. Note `test_rules_sphere.py:187` even carries the comment
`# S_8 reaches degree 5 in moment-condition rules; 4π/3 is degree 2.` — i.e. the
author *knew* the implemented rule is not a moment-condition rule and asserted
only degree 2, but the tag still claims degree 7.

Note also that the project already knows this in one place:
`.claude/rules/vv-testing.md` says "LS rules are NOT [exact] — they're optimised
for moment integration in transport, not arbitrary SH products". The rule file and
the code tag contradict each other.

### F2 `[M]` The `min()` in `product_mu_phi` is **numerically correct and sharp** — the *docstring's* incommensurable-units framing is the defect, not the arithmetic

Claimed `degree_of_exactness = min(2*n_mu - 1, n_phi - 1)`
(`rules_product.py:148`). Measured against the full monomial sweep on `S²`:

```
  n_mu= 1 n_phi=  1  claimed=  0  MEASURED=  0  OK
  n_mu= 1 n_phi=  4  claimed=  1  MEASURED=  1  OK
  n_mu= 2 n_phi=  2  claimed=  1  MEASURED=  1  OK
  n_mu= 2 n_phi=  4  claimed=  3  MEASURED=  3  OK
  n_mu= 2 n_phi=  8  claimed=  3  MEASURED=  3  OK      (polar-limited)
  n_mu= 2 n_phi= 16  claimed=  3  MEASURED=  3  OK      (polar-limited)
  n_mu= 3 n_phi=  4  claimed=  3  MEASURED=  3  OK
  n_mu= 3 n_phi=  6  claimed=  5  MEASURED=  5  OK
  n_mu= 4 n_phi=  4  claimed=  3  MEASURED=  3  OK      (azimuth-limited)
  n_mu= 8 n_phi=  4  claimed=  3  MEASURED=  3  OK      (azimuth-limited)
  n_mu= 4 n_phi=  8  claimed=  7  MEASURED=  7  OK
  n_mu= 4 n_phi= 16  claimed=  7  MEASURED=  7  OK
  n_mu= 6 n_phi= 12  claimed= 11  MEASURED= 11  OK
  n_mu= 8 n_phi= 16  claimed= 15  MEASURED= 15  OK
```

The tag is **exact and sharp on both sides** — the polar-limited and
azimuth-limited branches are each attained. The reason the `min()` is legitimate
is a *mapping the docstring never states*: for a monomial `x^a y^b z^c` restricted
to `S²`, the azimuthal integrand is `cos^aφ · sin^bφ`, whose **maximum
trigonometric frequency is exactly `a+b ≤ d`**. So "trig degree ≤ n_phi−1" and
"polynomial degree ≤ 2n_mu−1" are commensurable *after* the substitution
`trig-frequency = a+b`, and the composite bound is the min. The two units are
reconciled by the target function space (`Poly_d(R³)|_{S²}`), which is neither of
the two factor spaces.

⇒ The concrete defect is that **`degree_of_exactness` is an integer with no
attached subspace**, and the docstring papers over the reconciliation with the
word "conservative". Attaching a subspace would *document* this rather than break
it (see §6).

### F3 `[M]` `Chebyshev-Lobatto` in the CP kernel is a **naming falsity** — the nodes are Chebyshev-**Gauss** (first kind, interior)

`orpheus/derivations/continuous/flat_source_cp/geometry.py:96`
```python
_KI3_DEG: int = 63          # 64 Chebyshev-Lobatto nodes
```
`docs/theory/methods/collision_probability.rst:1965` — "at Chebyshev-Gauss-Lobatto
nodes". But the constructor is `np.polynomial.Chebyshev.interpolate` →
`np.polynomial.chebyshev.chebinterpolate`, whose numpy docstring reads
*"Interpolate a function at the Chebyshev points of the **first kind**"* — i.e.
`chebpts1`, which **excludes the endpoints** `[M]`:

```
  chebpts1(8) [Gauss, 1st kind]:   [-0.980785 ... 0.980785]     endpoints EXCLUDED
  chebpts2(8) [Lobatto, 2nd kind]: [-1.  -0.900969 ... 1.]      endpoints INCLUDED
```

The module's own docstring at `geometry.py:38` says "64 Chebyshev-Gauss nodes" —
correct — so the SAME FILE carries both the right and the wrong name 58 lines
apart. Three spellings for one node set (`Chebyshev-Gauss`,
`Chebyshev-Lobatto`, `Chebyshev-Gauss-Lobatto`) across code + theory page. This is
exactly the RANGE/SPACING/RULE conflation the campaign is chasing, in miniature:
Gauss-vs-Lobatto **is** the endpoint-inclusion axis, and it is mis-stated.

---

## 2. THE INVENTORY — every quadrature / node-set construction in the tree

Grouped by owning module. Column meanings: **Domain** = what is discretised;
**Nodes** = placement rule; **Weights**; **Exactness claim** = what the source
*says*, with the function space it says it on; **Adjudication** = my measurement
or reading.

### 2.1 `orpheus/numerics/` — the production quadrature layer

| # | `file:line` | Domain | Nodes | Weights | Exactness CLAIM (+ space) | Adjudication |
|---|---|---|---|---|---|---|
| N1 | `numerics/quadrature/rules_1d.py:79` | interval `μ∈[-1,1]` (polar cosine) | `np.polynomial.legendre.leggauss(n)` — Legendre roots | Christoffel numbers, `Σw=2` | `2n-1`, **polynomials** (Stoer–Bulirsch Thm 3.6.20, cited) | `[M]` **CORRECT and sharp** for n∈{1,2,3,4,6,8,16} |
| N2 | `numerics/quadrature/rules_product.py:111` | polar factor of `S²` | same `leggauss(n_mu)` | Christoffel | `2n_mu-1`, polynomials | `[M]` correct (inherits N1) |
| N3 | `numerics/quadrature/rules_product.py:115-116` | **circle `φ∈[0,2π)`** | `np.linspace(0, 2π, n_phi, endpoint=False)` — **left-endpoints** | **`2π/n_phi`, uniform** | `n_phi-1`, **trigonometric polynomials** | `[M]` correct: exact for `1≤|k|≤n_phi-1`, fails at `k=n_phi`. **THE FUSED EXPRESSION** — see §3–§5 |
| N4 | `numerics/quadrature/rules_product.py:148` | `S²` (the product) | outer `μ` × inner `φ` | `w_GL(μ)·2π/n_phi`, `Σw=4π` | `min(2n_mu-1, n_phi-1)`, "**conservative**" | `[M]` **CORRECT and sharp on both branches** — but the units reconciliation is undocumented (F2) |
| N5 | `numerics/quadrature/rules_sphere.py:93` | `S²` | `scipy.integrate.lebedev_rule(order)` — tabulated `O_h` orbits | tabulated, `Σw=4π` | `order`, polynomials (Lebedev 1976) | `[M]` **CORRECT and sharp** for order ∈{3,5,7,9,11,13,17} |
| N6 | `numerics/quadrature/rules_sphere.py:166-201` | `S²` | `μ²_p = μ²_1 + pΔ` recursion, then all 8 sign-octants of `(η,ξ,μ)` | **`4π/(8·n_octant)`, EQUAL within an octant** | `sn_order-1`, "conservatively", polynomials (cites Carlson–Lathrop 1968) | `[M]` **FALSE for N≥6 — actual degree is 3 for every N** (F1) |
| N7 | `numerics/measure.py:1022` | `[-1,1]` | `leggauss(n)` | Christoffel | `2n-1`, polynomials | `[M]` correct. **DUPLICATE of N1** — see §9 G1 |
| N8 | `numerics/measure.py:1074-1076` | `[-1,1]` with weight `(1-x²)^{-1/2}` | `cos((2i-1)π/2n)` — **Chebyshev–Gauss, 1st kind, closed form** | `π/n` uniform, `Σw=π` | `2n-1` **in the WEIGHTED sense** | `[M]` correct weighted; unweighted degree is **−1** (integrates nothing unweighted) — docstring says this explicitly and correctly |
| N9 | `numerics/measure.py:1126-1127` | arbitrary `[a,b]` | `a+(i+½)h` — **equispaced MIDPOINTS** | `h=(b-a)/n` uniform | `1` (midpoint rule) | `[M]` correct. **This is the unused "correct" sibling of N3** — see §4 |
| N10 | `numerics/roots_of_unity.py` **(untracked, in-flight)** | circle, `q`-th roots of unity | **group action**: integer quarter-split + sign flips + octant swap; trig evaluated on the first octant only | (node generator only — carries no weights) | none (it is a NODE generator, not a rule) | Reading only. Its docstring's measured claims are quoted in §4 |

### 2.2 `orpheus/moc/`

| # | `file:line` | Domain | Nodes | Weights | Exactness CLAIM | Adjudication |
|---|---|---|---|---|---|---|
| M1 | `moc/quadrature.py:89` | **half-circle `φ∈[0,π)`** | `np.linspace(0, π, n_azi, endpoint=False) + π/(2·n_azi)` — **equispaced MIDPOINTS** | `1/n_azi` uniform, **`Σω=1`** | **NONE STATED** anywhere in code or theory page | `[M]` It is the periodic trapezoid on the *quotient* circle `R/πZ`: exact for **π-periodic** trig polys of degree ≤ `n_azi-1`; only **algebraically** convergent for genuinely 2π-periodic content (residual 0.6533→0.6472→0.6376 → 2/π at n=4,5,16). See §3 |
| M2 | `moc/quadrature.py:20-33` `_TY_TABLES` | polar `θ` per half-space, `n_polar∈{1,2,3}` | **TABULATED literals** (`sin θ_p`) from Yamamoto 2007 Table 2 | tabulated, `Σ=0.5` per hemisphere | "optimised for the **Bickley-function integrals**" (module docstring line 4-6) | Reading only. **The only rule in the tree whose exactness space is NOT a polynomial/trig space** — it is a *minimax* fit on `Ki₃`. Theory page `method_of_characteristics.rst:350-354` restates it and adds "comparable to GL with 12–16 points" — no degree, no error bound in code |

### 2.3 `orpheus/cp/`

| # | `file:line` | Domain | Nodes | Weights | Exactness CLAIM | Adjudication |
|---|---|---|---|---|---|---|
| C1 | `cp/solver.py:170` | radial impact parameter `y ∈ [0, R]` | `composite_gauss_legendre(mesh.edges, n_quad_y)` — GL panels on the **mesh edges** | GL, per panel | inherits GL `2n-1` per panel | Reading only. `n_quad_y: int = 64` default (`cp/solver.py:66`). This is the *only* production consumer of the derivations quadrature contract |
| C2 | `cp/solver.py:94` | `μ ∈ [0,1]` (E₃ kernel) | none — closed form via `scipy.special.expn` | — | — | Not a quadrature (exact special function) |

### 2.4 `orpheus/mc/` — sampling, not quadrature

| # | `file:line` | Domain | Nodes | Weights | Adjudication |
|---|---|---|---|---|---|
| MC1 | `mc/solver.py:410-411` | direction | `theta = np.pi*rng.random()`, `phi = 2π*rng.random()` | implicit `1/N` | **Random, not a rule.** ⚠ Flag for a separate look: uniform-in-`θ` is **not** isotropic on `S²` (isotropic needs uniform in `μ=cos θ`); this is a 2-D solver using only `dir_x, dir_y = sinθ·cosφ, sinθ·sinφ`, so the sampled `sinθ` density is `∝1` rather than the correct `∝ sinθ/√(1-...)` — **outside this survey's scope, but it is a node-placement claim that does not hold.** Not filed (per brief) |
| MC2 | `mc/solver.py:289-290` | source position | `rng.random()*pitch` | — | Uniform-in-area sampling on a square: correct |

### 2.5 `orpheus/derivations/` — the semi-analytical verification machinery

This is where the *variety* lives. All of it returns the `Quadrature1D` contract
(`derivations/common/quadrature.py:66`) or is adaptive.

| # | `file:line` | Domain | Nodes | Weights | Exactness CLAIM (+ space) | Adjudication |
|---|---|---|---|---|---|---|
| D1 | `derivations/common/quadrature.py:249` `gauss_legendre(a,b,n)` | **arbitrary `[a,b]`** | affine-mapped `leggauss(n)`; `dps>53` routes to `mpmath.gauss_quadrature` | `½(b-a)·w` | `2n-1`, polynomials | Correct by construction (affine invariance) |
| D2 | `derivations/common/quadrature.py:268` `gauss_legendre_visibility_cone` | `[a,b]` with a **√-endpoint singularity** | `√(a²+u²Δ²)` for `u` = GL on `[0,1]` — a **Jacobian-absorbing substitution** | `u_w·u·Δ²/pts` | "**becomes spectral**" for integrands carrying `√(y²-a²)`; cites §22.7 with a Bernstein-ellipse analysis | Reading only. **A rule whose exactness space is a *transformed* polynomial space** — the strongest existing evidence that "the exactness space is a first-class parameter" |
| D3 | `derivations/common/quadrature.py:332` `composite_gauss_legendre` | sorted breakpoint list | GL panels | per panel | `2n-1` **per panel** | Correct; carries `panel_bounds`/`panel_sizes` so a consumer can index per-panel |
| D4 | `derivations/common/quadrature.py:500` `gauss_laguerre` | **`[0, ∞)`** with weight `e^{-x/σ}` | `np.polynomial.laguerre.laggauss(n)` (or mpmath) | `σ·w` | `2n-1` for `f` polynomial **after** the substitution `u=x/σ` | Reading only. Docstring: "*Currently used in diagnostics*" — a **shipped-but-unconsumed** production rule |
| D5 | `derivations/common/quadrature.py:381` `AdaptiveQuadrature1D` / `adaptive_mpmath` | `[a,b]` + breakpoint hints | **NO FIXED NODES** — `mpmath.quad` chooses at eval time | — | none (adaptive) | The type explicitly refuses to expose nodes; §"Notes" documents the static-vs-adaptive split |
| D6 | `derivations/common/quadrature_recipes.py:66` `chord_quadrature` | impact parameter `h∈[0,R]` | GL + **vis-cone panels at shell radii** | mixed | inherits D2/D3 | Geometry-aware recipe |
| D7 | `derivations/common/quadrature_recipes.py:188` `observer_angular_quadrature` | `ω ∈ [ω_low, ω_high]` | composite GL with **kink breakpoints at `arcsin(r_k/r_obs)` and their `π-ω` mirrors** | GL per panel | "recovers **spectral convergence on each smooth sub-panel**"; without it `O(Q^{-3/2})` | Reading only. **The range is a caller argument** — see §3 |
| D8 | `derivations/common/quadrature_recipes.py:305` `surface_centred_angular_quadrature` | `φ ∈ [φ_low, φ_high]` | composite GL + chord-quadratic kinks | GL per panel | as D7 | Range is a caller argument |
| D9 | `derivations/continuous/peierls_nystrom/geometry.py:419-428` `angular_range` | — | — | — | — | **`(-1,1)` for slab-polar, `(0,π)` for cylinder/sphere.** A *named, polymorphic* RANGE property on a geometry type. **The RANGE axis is ALREADY abstracted here** — see §3 |
| D10 | `derivations/common/kernels.py:262` `ki_n_mp` | `θ∈[0,π/2] → u∈[0,∞)` via `u=tanθ` | `mpmath.quad` adaptive on `[0, ∞)` | — | "removes the essential singularity at `θ=π/2`" | Substitution + adaptive |
| D11 | `derivations/common/kernels.py:307` | `u∈[0,∞)` | `scipy.integrate.quad` adaptive | — | `epsabs=1e-15, epsrel=1e-13` | **Second, independent realization of the SAME Bickley integral** as D10 |
| D12 | `derivations/continuous/flat_source_cp/geometry.py:96-118` | `τ∈[0,50]` | `np.polynomial.Chebyshev.interpolate(deg=63)` | (interpolant, not a quadrature) | comment: "*64 **Chebyshev-Lobatto** nodes*"; module docstring line 38: "*64 **Chebyshev-Gauss** nodes*"; theory page: "*Chebyshev-Gauss-**Lobatto***" | `[M]` **The nodes are Chebyshev-GAUSS (1st kind, endpoints EXCLUDED)** — F3. **Third realization of the Bickley integral** |
| D13 | `derivations/common/continuous_reference.py:323` | per-cell `[x_{i-½}, x_{i+½}]` | affine `leggauss(n_quad)` per cell | `½w·h` | "exact cell averages **to GL quadrature precision**" | Reading only. A 4th, hand-rolled composite-GL — see §9 G1 |
| D14 | `derivations/continuous/mms/moc.py:223-230` | 2-D cell `[0,P]²` | **tensor-product** `leggauss(32)⊗leggauss(32)` | `½P·w` each axis | none stated (`n_q = 32` hardcoded) | 5th hand-rolled GL |
| D15 | `derivations/discrete/sn/balance.py:367-368` | **circle `φ∈[0,2π)`** | `np.linspace(0, 2π, n_phi, endpoint=False)` | `2π/n_phi` | none (it is a *verification* of `Σwη=0`) | **VERBATIM TWIN of N3** re-derived inside a derivation module — see §9 G1 |
| D16 | `derivations/discrete/sn/dsa.py:123-158` | slab angular | **symbolic** S4 quadrature; moments stay symbolic | symbolic | "quadrature-**INDEPENDENT**" (line 254) | The one place the *absence* of a concrete rule is the point |
| D17 | `derivations/continuous/singular_eigenfunction/core/half_range.py:99` | `t∈[0,1]` | `np.linspace(0,1,n_grid)` — **equispaced, endpoint-INCLUSIVE** | (a grid for an X-function evaluation, not a weighted rule) | none | A 3rd `linspace` convention (`endpoint=True`) |
| D18 | `derivations/continuous/fn_method/moment_space.py:508,557` | `z∈[-a,a]` / `r∈[0,R]` | `np.linspace` evaluation panels | — | none | Evaluation grids, not rules |

### 2.6 `orpheus/geometry/`, `diffusion/`, `homogeneous/`, `transport/`, `fuel/`, `thermal_hydraulics/`, `kinetics/`

- **`geometry/mesh.py`, `geometry/factories.py`** — `np.linspace` for **spatial** cell edges. Node placement, but on the *spatial* axis, and the mesh already owns that concept as a first-class type (`Mesh1D.edges`). Out of the angular-quadrature scope; listed for completeness.
- **`diffusion/`** — `[M]` no angular quadrature at all (grep for `linspace|leggauss|quad` returns nothing in `orpheus/diffusion/`). Correct: P1/diffusion has integrated the angle away analytically.
- **`homogeneous/`** — `[M]` likewise none. 0-D infinite medium.
- **`transport/`** — consumes `Quadrature`; constructs none. (`transport/operators/scattering.py:483` reads `quadrature.spherical_harmonics(...)`, `:514` `quadrature.angular_frame(L)`.)
- **`fuel/solver.py:180,193`, `thermal_hydraulics/`, `kinetics/`** — `np.linspace` radial nodes + `scipy.integrate.solve_ivp`. ODE integration, not quadrature.

### 2.7 Usage census `[M]`

```
   Quadrature.gauss_legendre   orpheus/: 13   tests/: 450
   Quadrature.lebedev          orpheus/:  4   tests/: 137
   Quadrature.level_symmetric  orpheus/:  4   tests/: 196
   Quadrature.product          orpheus/:  6   tests/: 117
```

The four families are all live; `level_symmetric` (the one with the false
exactness tag) has 196 test call sites.

---

## 3. THE RANGE AXIS — every range choice, stated or implicit

| Range | Where | Stated or implicit? |
|---|---|---|
| **`[0, 2π)` full circle** | N3 `rules_product.py:115`; D15 `balance.py:367` | **Implicit.** The code comment says only "*matching the legacy ProductQuadrature contract pinned by the regression snapshots*" — a *provenance* justification, not a mathematical one |
| **`[0, π)` half circle** | M1 `moc/quadrature.py:89` | **STATED, and justified.** `moc/quadrature.py:44-47`: "*the supplementary angle (phi + pi) is the same physical track traversed in the opposite direction*"; theory page `method_of_characteristics.rst:328-332` repeats it. **This is the one place a symmetry quotient is explicitly taken and explicitly argued.** But the *quotient* is never named as such: there is no `HalfRange`, no `quotient_by`, no π-periodicity type. The consequence — that the rule is spectral only on π-**periodic** integrands — is nowhere stated `[M]` |
| **`[-1, 1]` interval (polar cosine)** | N1, N2, N7, N8, N9, D9(slab) | Stated in every docstring |
| **`[0, 1]` half interval** | D17 `half_range.py:99`; `peierls_nystrom/geometry.py:2464,2787` (`gauss_legendre(0.0, 1.0, n_quad)`) | Implicit — a bare literal pair |
| **`[0, π]` (`ω`, curvilinear ray angle)** | D9 `peierls_nystrom/geometry.py:428` | **STATED AND NAMED.** `angular_range` is a `@property` on the geometry type that *returns* `(-1.0, 1.0)` for slab-polar and `(0.0, np.pi)` otherwise. **The RANGE is already a first-class, polymorphic concept — in the derivations tree only** |
| **arbitrary `[a, b]`** | D1, D3, D6, D7, D8, N9 | A **caller argument** in every case |
| **`[0, ∞)`** | D4 `gauss_laguerre`; D10, D11 (Bickley) | Stated |
| **`S²` full sphere** | N4, N5, N6 | Stated |
| **hemisphere / octant** | not a range anywhere — it is a `restrict()`/`partition_by()` **operation** on a full-sphere measure (`measure.py:597,652,750`, `half_range_clean` flag) | The half-sphere concept exists only as a *post-hoc filter*, never as a rule's domain |

**What the axis shows.** The RANGE is *already* a named, polymorphic property in
exactly one place (D9, on a derivations geometry type) and a hard-coded literal
everywhere else. MoC's `[0, π)` is the sharpest evidence that a range choice is a
*modelling* decision: it is a **symmetry quotient** (`R/πZ`) taken for a stated
physical reason, and it changes the rule's exactness space from
"2π-periodic trig polys" to "π-periodic trig polys" — a fact `[M]` measured
above and stated nowhere.

**And the range choice is currently a live BUG.** Issue **#326** (`module:sn`,
`module:geometry`, `level:L1`, `type:bug`, open) says exactly this about N3:

> "The key is `η = sinθ·cos φ`, which is **2-to-1 over φ ∈ [0, 2π)**: the azimuthal
> mirror pair `(φ, 2π−φ)` shares `η` and differs only in `ξ = μ_y`. So the level was
> never totally ordered by this key."

with measured consequence 0.4–12 % on heterogeneous cylindrical problems. The
cylindrical *level* is a half range that the product rule realises as a full
circle.

**A theory page currently argues the opposite.**
`docs/theory/methods/sn/curvilinear_one_group.rst:3891-3956`
(`sn-direct-seed-circle-vs-interval`) makes the full circle load-bearing:

> "**Cylinder — the redistribution axis is a circle.** … lives on a **circle**
> `[0, 2π)` — a *periodic* domain. … Equispaced sampling of a smooth periodic
> function is the trapezoidal rule, which on a circle is **spectrally accurate** …
> there is no accuracy penalty for the choice."
> "The principle in one line: **a periodic redistribution axis gives edge-inclusion
> for free; an interval axis makes you pay for it with a separate seed.**"

That page's premise is the very thing #326 calls a double cover. **Whichever way
the #326 ruling lands, this section needs rewriting** — it is the durable
statement of the RANGE choice in the corpus and it is currently in tension with
an open L1 bug. (Consistent with my prior note that this section's periodicity
principle is falsified by the half-range analysis.)

---

## 4. THE SPACING AXIS — equispaced vs Gauss vs tabulated

### Where equispaced is LOAD-BEARING (equal weights are optimal, not a compromise)

- **N3 (SN azimuth on `[0,2π)`).** Load-bearing *given the range*. `[M]` The
  periodic trapezoid is exact for `1≤|k|≤n_phi−1` and, for smooth periodic
  integrands, converges faster than any power of `1/n` (the theory page states
  this at `curvilinear_one_group.rst:3898-3900`). Equal weights here are the
  **Gauss rule of the circle** — no other node placement improves on it for a
  periodic integrand. Genuinely optimal.
- **M1 (MoC azimuth on `[0,π)`).** Load-bearing *only if the integrand is
  π-periodic*, which it is **after** the forward/backward track pairing
  (`moc/core.py:128-136` sums both sweep directions into one `omega_azi` weight).
  `[M]` measured: exact for π-periodic freq ≤ n_azi−1; residual for a 2π-periodic
  `e^{iφ}` is `6.53e-1 → 6.38e-1` at n=4→16, i.e. only algebraic. The optimality
  is *conditional on the quotient*, and nothing in the code records that condition.

### Where equispaced is INCIDENTAL

- **N9 `equispaced(a,b,n)`** (`measure.py:1085`) — a general midpoint rule with
  `degree_of_exactness=1`. Not periodic; equal weights are just the midpoint rule.
  **`[M]` It has ZERO consumers in `orpheus/`** — grep for `equispaced(` outside
  `measure.py` returns nothing in production. Its own docstring (`measure.py:1100-1107`)
  says it exists as the *would-be* azimuthal primitive and that "*the project's
  existing code uses left-endpoints, but this primitive offers midpoints because
  they integrate constants exactly while preserving symmetry under reflection
  through the centre of the interval*". **A shipped, documented, unconsumed
  alternative realization of the same axis.**
- **D17 `linspace(0,1,n_grid)`** with `endpoint=True` — an evaluation grid.
- **D18, and every `geometry/` spatial `linspace`** — evaluation/mesh grids.

### The OFFSET is a third sub-axis, and it is invisible

Three different offsets of the same equispaced family ship today:

| Offset | Where | Why |
|---|---|---|
| left-endpoint (`c=0`) | N3, D15 | "matching the legacy contract" (provenance) |
| midpoint (`c=½`) | M1, N9 | M1: unstated. N9: "*integrate constants exactly while preserving symmetry under reflection through the centre*" |
| endpoint-inclusive | D17 | unstated |

`[M]` **The offset does NOT change exactness on the circle:**

```
    n= 4 left-endpoint (SN product)   exact for 1<=|k|<= 3  (claim n-1=3)
    n= 4 midpoint                     exact for 1<=|k|<= 3  (claim n-1=3)
    n= 5 left-endpoint (SN product)   exact for 1<=|k|<= 4  (claim n-1=4)
    n= 5 midpoint                     exact for 1<=|k|<= 4  (claim n-1=4)
    n= 8 left-endpoint (SN product)   exact for 1<=|k|<= 7  (claim n-1=7)
    n= 8 midpoint                     exact for 1<=|k|<= 7  (claim n-1=7)
```

So the offset is a pure **symmetry / degeneracy-avoidance** choice, not an
accuracy one — and it is currently made three different ways with three
different (or absent) justifications. This is the sharpest single argument that
the three fused choices in `rules_product.py:115-116` are genuinely independent:
changing the offset changes *which symmetries the node set exactly respects*
while leaving the exactness tag identical.

### The EXACT-GENERATION sub-axis (in flight)

`orpheus/numerics/roots_of_unity.py` *(untracked, in-flight; issue **#325**,
`level:L0`)* is a fourth realization of "place `n` equispaced nodes on a circle" —
same nodes to `≤ 9.99e-16`, but generated by the **group action** instead of by
evaluating `cos`/`sin`. Its docstring's measured claims (read, not re-measured
by me):

> "x-mirror `k -> (n-k) mod n` bit-exact: this module: every q / linspace+cos: NO q"
> "y-mirror `k -> (n/2 - k)` bit-exact: this module: every even q / linspace+cos: NO q"
> "45-diagonal `|cos| == |sin|` bit-exact: this module: every `q % 8 == 0` / linspace+cos: NO such q"

Issue #325's own table separates the families exactly along this axis:
`lebedev` and `level_symmetric` are **algebraic** (sign/permutation orbits) and
exact; `product` and **MoC azimuthal** evaluate a parameterization and are not.

`[M]` I confirmed the consequence is observable through the registry's own flags:
for **odd `n_mu`** the product rule *should* place 4 nodes on coordinate axes
(GL has a node at `μ=0`, and `n_phi%4==0` hits `φ∈{0,π/2,π,3π/2}`), but only
**1** of the 4 has exact zeros:

```
  product(3,4): on-axis EXACT=1  on-axis tol1e-14=4  N=12
  product(3,8): on-axis EXACT=1  on-axis tol1e-14=4  N=24
  product(5,4): on-axis EXACT=1  on-axis tol1e-14=4  N=20
```

(The 1 that is exact is `φ=0`, where `sin(0)==0.0` exactly — precisely the
asymmetry #325 describes.)

---

## 5. THE RULE AXIS — the actual rule and its exactness space, per 1-D factor

| 1-D factor | Rule | Exactness SPACE | Realizations in-tree |
|---|---|---|---|
| polar `μ∈[-1,1]` | **Gauss–Legendre** (open) | `P_{2n-1}` on `[-1,1]`, Lebesgue | N1, N2, N7, D1, D3, D13, D14 — **7 spellings, one rule** |
| polar `μ` | **Gauss–Lobatto** (closed, endpoints included) | `P_{2n-3}` | **ZERO in-tree.** Studied and **declined** — see below |
| polar `μ` (weighted) | **Gauss–Chebyshev 1st kind** | `P_{2n-1}` w.r.t. `(1-x²)^{-1/2}` | N8 — shipped, **unconsumed in `orpheus/`** |
| `[0,∞)` weighted | **Gauss–Laguerre** | `P_{2n-1}` after `u=x/σ`, w.r.t. `e^{-u}` | D4 — shipped, "*currently used in diagnostics*" |
| azimuth `φ` (circle) | **periodic trapezoid** (= equispaced+equal weights) | trig polys of degree `≤ n-1`; **spectral** on smooth periodic | N3, M1, D15, N9(midpoint variant) |
| azimuth `φ` (circle) | Gauss / Chebyshev / Lobatto / Clenshaw–Curtis / Simpson | — | **ZERO.** No second rule on the circle exists or is proposed anywhere |
| `[a,b]` with a `√` endpoint singularity | **substituted GL** (Jacobian-absorbing) | polynomials in the *transformed* variable `u` | D2 — **the only rule whose exactness space is explicitly a transformed one** |
| polar `θ` for `Ki_n` | **Tabuchi–Yamamoto** (tabulated) | **minimax on the Bickley family** `Ki₃`, not a polynomial space | M2 |
| `τ` interpolation | **Chebyshev interpolation** (1st-kind nodes) | near-minimax `P_63` on `[0,50]` | D12 |
| adaptive | `mpmath.quad` (tanh-sinh/GL adaptive) | none — error-controlled | D5, D10 |
| adaptive | `scipy.integrate.quad` (QUADPACK) | none — error-controlled | D11 |
| **QMC** | Owen-scrambled Sobol′ | Koksma–Hlawka on **bounded Hardy–Krause variation** | **ZERO in production**; validated in a diagnostic, issue **#128** |

### Is there a SECOND realization of "a rule on the circle"? — NO, and that is the finding

Every circle rule in the tree is the same periodic trapezoid, differing only in
range (`2π` vs `π`) and offset (`0` vs `½`). There is no Gauss-on-the-circle, no
Clenshaw–Curtis, no Simpson, and nothing in the issue tracker proposes one.
**On the mathematics, this is correct**: for a smooth periodic integrand the
periodic trapezoid *is* the optimal rule (it is spectrally accurate and its
"Gauss" analogue on the circle is itself). A second *rule* on the circle would be
a pessimisation.

⇒ *Inference:* the circle's RULE slot is genuinely **single-realization** and
should stay a property, not a type. The variation on the circle is entirely in
**range** and **offset/generation** — which is where ≥2 realizations DO exist.

### The Gauss–Lobatto question — LIVE as a documented ruling, ARCHAEOLOGY as code

`docs/theory/methods/sn/curvilinear_one_group.rst:3958-4026`
(`sn-direct-seed-lobatto-study`). Verbatim, the operative sentences:

> "**Affordable.** At resolved angular order (`N ≥ 8` …) Gauss–Lobatto tracks
> Gauss–Legendre at a bounded `~1.2×` error penalty … fine-`N` GL and GLob agree to
> `< 6` pcm — the two rules converge to the **same** `N → ∞` transport limit"
> "**But not a drop-in.** A pole node lands on the level's lower edge, so the
> first-ordinate weight `τ_{raw,0} = 0` — and the production Morel–Montry
> recurrence **divides its first step by that weight** … separately, the R12a
> presence predicate keys on `τ_raw ∈ (0,1)`"
> "**Ruling: affordable but architecturally declined.** … the reason is
> architectural, not numerical."

And the note:

> "The Gauss–Lobatto study is a set of scratch diagnostics
> (`scratch/experimental/glob_sphere_study/` and
> `derivations/diagnostics/diag_glob_0{1..5}_*.py` — 33 green diagnostics …).
> They are **uncommitted** and are promotion targets *only if* a pole-node scheme
> is ever adopted; do not promote them otherwise."

**Reading for this survey:** Gauss–Lobatto is a *fully-characterised, measured,
declined* second realization of the interval rule. The decline reason is **not**
"we don't need a second rule" — it is that adopting it forces a downstream
architectural change (`τ_raw=0` breaks a recurrence and a structural predicate).
That is precisely the shape of "the RULE choice is not encapsulated": swapping
the rule leaks into the solver. If the RULE were a first-class abstraction with a
declared `includes_endpoints` property, the recurrence and the R12a predicate
would read that property instead of a magic `τ_raw ∈ (0,1)` interval test.

---

## 6. THE `degree_of_exactness` AUDIT — producers, propagators, consumers

### 6.1 Producers (8)

| Producer | Value | Space it *means* | Verdict |
|---|---|---|---|
| `quadrature/rules_1d.py:85` | `2n-1` | polynomials on `[-1,1]` (Lebesgue) | `[M]` true |
| `quadrature/rules_product.py:148` | `min(2n_mu-1, n_phi-1)` | polynomials on `S²` — via an unstated `trig-freq = a+b` mapping | `[M]` true & sharp; framing wrong (F2) |
| `quadrature/rules_sphere.py:102` | `order` | polynomials on `S²` | `[M]` true |
| `quadrature/rules_sphere.py:278` | `sn_order-1` | polynomials on `S²` | `[M]` **FALSE for N≥6** (F1) |
| `measure.py:1027` | `2n-1` | polynomials on `[-1,1]` | `[M]` true. **Duplicate of rules_1d; `[M]` zero production consumers** (only `tests/numerics/test_measure.py`) |
| `measure.py:1081` | `2n-1` | polynomials **w.r.t. the weight `(1-x²)^{-1/2}`** — a DIFFERENT integral | `[M]` true weighted, **`-1` unweighted**. Same integer, incompatible meaning |
| `measure.py:1132` | `1` | polynomials on `[a,b]` | `[M]` true |
| `tests/sn/_test_helpers.py:1039` | `min(2*n_mu-1, n_phi-1)` | — | **A test-side twin producer** re-spelling the product formula |

Note the third column: **four different function spaces already share one
integer field** (Lebesgue-`[-1,1]` polys / `S²` polys / Chebyshev-weighted polys /
`[a,b]` polys). `measure.py:1081` and `measure.py:1027` both say `2n-1` and mean
integrals that differ by the weight function. Nothing in the type distinguishes
them.

### 6.2 Propagators (the `DiscreteMeasure` composition algebra)

| `file:line` | Operation | Rule |
|---|---|---|
| `measure.py:467-470` | `__mul__` (tensor product) | `min(p_μ, p_ν)` if both known, else `None` |
| `measure.py:528-531` | (2nd composition op) | `min(p_μ, p_ν)` |
| `measure.py:597,617` | `pushforward` | **dropped to `None`** unless `φ` is the identity |
| `measure.py:652,667` | `restrict` | **dropped to `None`** |
| `measure.py:750,786` | `partition_by` | **dropped to `None`** |
| `measure.py:803-823` | `with_tags` | explicit caller override |

The propagator is `min()` — the same incommensurability as F2, one level up, and
here it is *not* rescued by a mapping: `gauss_legendre(5) * gauss_chebyshev(4)`
would report `min(9, 7) = 7` across a Lebesgue axis and a Chebyshev-weighted axis.
`tests/numerics/test_measure.py:219` pins `p.degree_of_exactness == 5` for such a
product.

The drop-to-`None` behaviour on `restrict`/`partition_by`/`pushforward` is
**correct and conservative**, and is the one place the current design admits it
does not know the answer.

### 6.3 Consumers — who READS the attribute and COMPUTES with it

`[M]` Exactly **one** site in the whole tree computes with the value:

```python
# tests/numerics/test_spherical_harmonic_basis.py:120-121
deg = measure.degree_of_exactness
L = deg // 2  # safe order so ℓ+ℓ' ≤ deg
```

It interprets the integer as **spherical-harmonic** exactness and halves it. That
folklore is documented once, in
`docs/theory/foundations/spherical_harmonics.rst:293`:

> "For a discrete angular cubature … whose ``degree_of_exactness`` is at least
> :math:`2L`, the no-prefactor real `Y_ℓ^m` satisfy the **discrete** orthogonality"

Every other `degree_of_exactness` reference in `tests/` is an **assertion pinning
the tag's value** (`test_rules_1d.py:43`, `test_rules_product.py:40`,
`test_rules_sphere.py:39,111`, `test_registry.py:190,207`,
`test_measure.py:55,87,219,315,318,349`, `test_measure_partition.py:237,243`) —
label pins, not property checks.

The **selector** (`registry.py:686`) does not read the attribute: it calls
`spec.degree_of_exactness_for(target)` — the *inverted formula*, hard-coded per
family in `_gl1d_invert` / `_lebedev_invert` / `_ls_sn_invert` / `_product_invert`.
So the formula is duplicated: once in the rule (as the tag) and once in the
registry (as the inverse), with no consistency check between them.

### 6.4 `[M]` The measured cost of the missing subspace — a production consequence

The discrete SH Gram (the object the P_L scattering source is built on) measured
directly:

```
  level_symmetric(8)   L=1  tag=  7  needs deg>= 2  max|Gram err| = 2.487e-14
  level_symmetric(8)   L=2  tag=  7  needs deg>= 4  max|Gram err| = 2.435e-01
  level_symmetric(8)   L=3  tag=  7  needs deg>= 6  max|Gram err| = 2.705e-01
  level_symmetric(16)  L=1  tag= 15  needs deg>= 2  max|Gram err| = 7.994e-14
  level_symmetric(16)  L=2  tag= 15  needs deg>= 4  max|Gram err| = 4.524e-01
  lebedev(11)          L=2  tag= 11  needs deg>= 4  max|Gram err| = 1.066e-14
  product(6,12)        L=2  tag= 11  needs deg>= 4  max|Gram err| = 5.329e-15
```

At `L=2` the level-symmetric discrete SH basis is **not orthogonal at the 25–45 %
level**, while its tag certifies degree 7 / 15. And the production frame does not
notice: `[M]` `quadrature.angular_frame(L).gram` is a
`SphericalHarmonicSpace` whose `metric_per_ell` is the **analytic** Gram —

```
  metric_per_ell (the frame ASSUMES this Gram): [12.56637061  4.1887902   2.51327412]
  analytic 4pi/(2l+1)                        : [12.566370614, 4.1887902048, 2.5132741229]
```

— i.e. the Galerkin frame assumes an orthogonality the quadrature does not
provide, and the only thing that was supposed to guarantee it is the false tag.
(`transport/operators/scattering.py:514` is the production call site.)

### 6.5 Would attaching a SUBSPACE break or fix the consumers?

| Consumer | Effect of attaching a subspace |
|---|---|
| `test_spherical_harmonic_basis.py:120` `L = deg // 2` | **FIXES.** The `//2` is folklore converting "poly degree on `S²`" to "SH degree"; with a subspace it becomes `measure.is_exact_on(SphericalHarmonics(L))` and stops being a magic halving |
| `registry` V-filter | **FIXES.** The filter's real question is "is this rule exact on the space my consumer integrates?" — currently spelled as an integer that means four things |
| `measure.__mul__` `min()` | **FIXES.** The honest composite space is the tensor product `V_μ ⊗ V_ν`; `min()` is the shadow it casts on a common scale. F2's reconciliation becomes a derived query rather than an unstated assumption |
| `rules_sphere.level_symmetric_sn` | **FIXES** — `sn_order-1` on `Poly(S²)` becomes unspellable; the honest tag is the (small) moment set it actually satisfies, or `None` |
| `measure.gauss_chebyshev` | **FIXES** — the "weight function is in the integral, not the quadrature" caveat (`measure.py:1052-1057`) becomes machine-readable instead of prose |
| **`moc.MOCQuadrature`** (M2, TY) | **ENABLES.** Today TY carries **no exactness field at all** because its space is a minimax fit on `Ki₃`, not a polynomial space. It is the one rule that *cannot* join the tagged registry under an integer field |
| **All the tag-pinning asserts** | **BREAK cosmetically** (they compare an int) — but they were pinning a label, and F1 is exactly the failure mode that survives label-pinning |

⇒ Attaching a subspace **breaks nothing that was checking a property** and fixes
every site that was doing arithmetic on the integer. The one *hard* case (TY) is
hard today too — it just fails silently by being absent.

---

## 7. FORESEEABLE NEEDS — documented vs inferred

### 7.1 DOCUMENTED (open GitHub issues, verbatim)

| Issue | Axis it lands on | The documented ask |
|---|---|---|
| **#325** `module:sn,moc,geometry` `level:L0` | **SPACING / generation** | "*generate symmetry-defined node sets from the group action, not by evaluating a parameterization*". Its own table splits the four families exactly on this axis: `lebedev`/`level_symmetric` **algebraic** (exact mirrors); `product`/**MoC azimuthal** parameterization-evaluated (no exact mirror for any `n`). Consequences it names: `TANGENTIAL_EPS = 4·eps` "*an epsilon deciding which ordinates are physically grazing — a magic number deciding physics*"; `reflection_index` degraded to a **nearest-neighbour search** |
| **#326** `level:L1` `type:bug` | **RANGE** | "*the azimuthal mirror pair `(φ, 2π−φ)` shares `η` … the level was never totally ordered by this key*" — measured **0.4–12 %** flux error on heterogeneous cylinders, discriminated by spatial heterogeneity, **not converging in `n_phi`** |
| **#128** `module:cp` | **RULE (a genuinely new family)** | Owen-scrambled **Sobol′ QMC** on the angular dimension. Validated in a diagnostic: "*32 scrambles × N=4096 … 95 % bootstrap CI widths … 20-100× tighter*". Proposes a literal `quadrature={"product_gauss","qmc_angular","qmc_3d"}` parameter. Exactness space: **bounded Hardy-Krause variation** (Koksma-Hlawka) — a fourth kind of space |
| **#109** `module:cp` | **RANGE + WEIGHT (change of variable)** | τ-coordinate: "*Gauss-Laguerre on `[0, τ_max]` with `e^{-τ}` weight is the ideal quadrature; expect 2-4× fewer nodes*". Also proposes §22.4 exp-stretched `μ = e^{-v}` and §22.5 **Gauss-Jacobi endpoint weights**. `gauss_laguerre` (D4) already ships for this and is unconsumed |
| **#123** `module:cp,tests` `level:L1` | **the ≥2-realization criterion itself** | Mandates "*≥ 2-quadrature signed-error stability protocol*": a closure claim is only verified if it holds "*at every quadrature in `quads`*", with no sign flip and monotone convergence. **This is the project's own, already-ratified statement that ≥2 realizations is what makes a claim provable** |
| **#145** | tooling | per-reference quadrature floor sweep + at-floor pytest marker |
| **#265** `module:data,geometry,numerics` | **the abstraction itself** | "Data-layer invariant value-objects (S13): **Quadrature laws**, albedo α∈[0,1], …" — a filed intent to make quadrature properties invariant-bearing value objects |
| **#235** `module:sn` `level:L1` | **RANGE / closure** | "*cylindrical anisotropic accuracy needs a 2-D (η,φ) angular closure*" |
| **#191** | consumption | retire `for n in range(quad.N)` for an `Ordinate` iterator |
| **#154 / #153 / #152 / #151** | — | the registry / `AngularQuadrature`→`DiscreteMeasure` / subgroup-lattice / `DiscreteMeasure` issues are **still OPEN although the code exists**. Hygiene note, not a need |

### 7.2 DOCUMENTED in code/docs (deferrals with the seam already cut)

- **Hexagonal / `D_6h`.** `numerics/symmetry.py:155`: "*``Dnh(6)`` for hexagonal
  lattices (**forward-looking — not yet consumed**)*". The **group is already
  implemented** (`symmetry.py:263`); what is missing is a rule invariant under it
  and a `GEOMETRY_GROUPS` entry. `registry.py:492`: "*New geometries (hexagonal
  lattice, 2-D / 3-D spherical, …) are added here, not in the selector itself.*"
  ⇒ A `D_6h`-invariant angular rule has **no in-tree realization** and no
  candidate: neither Lebedev (`O_h`) nor level-symmetric (`O_h`) nor product
  (`SO(2)`/`C_{n_φ}`) is `D_6h`-invariant in general. *Inference:* the natural
  construction is exactly the product rule with `n_φ` a multiple of 6 — i.e.
  **a RANGE/SPACING choice, resolved by making `n_φ` divisible by the lattice
  order**, which the current fused expression cannot express or check.
- **2-D / 3-D spherical SN.** `registry.py:119-124`: "*When 2-D / 3-D spherical SN
  lands, this table will gain a ``sphere2d`` / ``sphere3d`` entry tagged `O_h` or
  `O(3)`.*"
- **Gauss–Lobatto.** Characterised, measured, **declined** for architectural
  reasons (§5). The 33 diagnostics are uncommitted and explicitly *not* to be
  promoted unless a pole-node scheme is adopted.
- **Anisotropic SH requirement.** `docs/theory/foundations/spherical_harmonics.rst:293`
  states the `degree_of_exactness ≥ 2L` requirement. `.claude/rules/vv-testing.md`
  states which factories meet it: "*Lebedev order ≥ 2L, `product_mu_phi(n_μ≥L+1,
  n_φ≥2L+1)` are exact; **LS rules are NOT***".

### 7.3 INFERRED (my reading, marked as inference)

- *Inference:* **anisotropic scattering above `P_1` is currently unsound on the
  level-symmetric family.** Grounded in the `[M]` Gram measurement (§6.4): the
  discrete SH basis on `level_symmetric` is non-orthogonal at 25–45 % for `L ≥ 2`,
  and the frame assumes the analytic Gram. Whether any production configuration
  actually pairs them I did **not** verify at call-site precision — grep shows
  `Quadrature.level_symmetric` and anisotropy tokens co-occurring in ~20 test
  modules, but co-occurrence in a file is not a pairing in a call. **This is the
  highest-value thing to check next and I flag it as unverified.**
- *Inference:* the **`D_6h` need converts the azimuthal `n_φ` from a resolution
  knob into a group-order constraint** (`lattice_order | n_φ`). Today nothing can
  express or check that, because `n_φ` is consumed inside one `linspace`.
- *Inference:* #109's τ-coordinate and #128's QMC both need a rule to declare
  **which weight function / which variation class** it is exact on — i.e. both
  independently demand the same missing subspace field as §6.
- *Inference:* MoC and SN will eventually need to share one azimuthal primitive
  (the two differ only in range and offset). #325 already treats them as one
  problem ("MoC azimuthal" is a row in its table); #279/#261 push CP/MoC/SN onto
  shared operators generally.

---

## 8. SCORECARD — which choice has ≥2 genuine realizations?

Using the user's criterion: *a choice earns an abstraction when there are ≥2
genuine realizations in the codebase's real and foreseeable needs, because ≥2
realizations let a mathematical invariant be tested across them.*

| Axis | In-tree realizations | Foreseeable (documented) | Verdict |
|---|---|---|---|
| **RANGE** | `[-1,1]`, `[0,1]`, `[a,b]`, `[0,2π)`, **`[0,π)`**, `[0,π]`, `[0,∞)`, `S²` — **8**, and D9 already implements it as a polymorphic `angular_range` property | #326 needs the half range on the SN azimuth; #109 needs `[0,τ_max]` | **EARNS IT — strongest case.** ≥2 realizations exist *for the same physical angle* (SN's `[0,2π)` vs MoC's `[0,π)`), the choice is a **symmetry quotient**, and the invariant it enables is testable: *"the quotient rule and the full rule agree on the G-invariant sector and only there"* — precisely the #326 adjudicator |
| **SPACING (node placement)** | Gauss (7 spellings), Chebyshev-Gauss (2), Laguerre, equispaced (3 offsets), tabulated TY, adaptive (2 engines) — **easily ≥2** | QMC/Sobol′ (#128); Gauss-Jacobi (#109 §22.5) | **EARNS IT.** But note the sub-split: *within the circle*, only equispaced exists and should — the variation is in **offset** and **generation method** (#325), not in placement family |
| **RULE / weights on the CIRCLE** | **ONE** — the periodic trapezoid, everywhere | none proposed anywhere | **DOES NOT EARN A TYPE.** Mathematically it *should not* — the periodic trapezoid is the circle's Gauss rule. Keep it a property |
| **RULE / weights on an INTERVAL** | GL, Gauss-Chebyshev, Gauss-Laguerre, midpoint, adaptive — **5** | Gauss-Lobatto (studied, declined, fully characterised); Gauss-Jacobi (#109) | **EARNS IT**, and the Lobatto ruling shows the *cost of not having it*: the swap leaked into the Morel-Montry recurrence and the R12a predicate |
| **EXACTNESS SPACE** | `Poly[-1,1]`, `Poly(S²)`, `Poly` w.r.t. `(1-x²)^{-1/2}`, `Poly` w.r.t. `e^{-x}`, trig polys, π-periodic trig polys, transformed-variable polys (D2), Bickley-minimax (TY), *nothing* (adaptive) — **9 distinct spaces already**, all collapsed into one unlabelled `int` | Hardy-Krause variation (#128) | **EARNS IT — and it is the axis with the most in-tree realizations already.** It is also the one whose absence produced F1 (a 12-degree lie with a live consumer) |
| **NODE GENERATION (exact vs evaluated)** | 2 algebraic (`lebedev`, `level_symmetric`), 2 evaluated (`product`, MoC), 1 in-flight exact generator (`roots_of_unity`, untracked) | #325 | **EARNS IT** — a property/protocol on the node set, with the invariant "`node_set` is closed under `G` **bit-exactly**" testable across all four families (`SubgroupOfO3.is_invariant` is the existing checker; `roots_of_unity` is the generator) |

### The one-line reading

The fused expression `np.linspace(0, 2π, n_phi, endpoint=False)` + `2π/n_phi`
fuses **four** things, not three: RANGE, SPACING, RULE — **and the generation
method**. Of the four, three have ≥2 genuine realizations in-tree today
(range, spacing/offset, generation) and the fourth (the rule on the circle)
correctly has one and should stay a property. Sitting under all of them, the
axis with the *most* existing realizations and the only one with a **measured
falsehood shipping today** is the **exactness space**.

---

## 9. GAPS — expected but not found

- **G1 — Twin paths.** The GL rule is spelled **7 times** (N1, N2, N7, D1, D3,
  D13, D14) and the periodic-trapezoid circle rule **twice verbatim**
  (`rules_product.py:115-116` ≡ `derivations/discrete/sn/balance.py:367-368`).
  `measure.gauss_legendre` (N7) `[M]` has **zero production consumers** — it is a
  test-only duplicate of `rules_1d.gauss_legendre_on_mu` differing only in not
  setting `invariance_group`.
- **G2 — Shipped, documented, unconsumed primitives** `[M]`:
  `measure.equispaced` (0 production consumers), `measure.gauss_chebyshev` (0),
  `derivations.gauss_laguerre` (0, docstring says "*currently used in
  diagnostics*"). These are **already the second realizations** the criterion
  asks for — they exist but are not wired, so no invariant is tested across them.
- **G3 — MoC is outside the tagged system entirely.** `MOCQuadrature`
  (`moc/quadrature.py:36`) is a plain frozen dataclass with **no
  `DiscreteMeasure`, no `invariance_group`, no `degree_of_exactness`, no
  registry entry**. Its weights sum to **1** (a probability normalisation) while
  every `numerics` sphere rule sums to `4π` — two normalisation conventions with
  no crosswalk. Its TY table is the only exactness claim in the tree that cannot
  be expressed as a polynomial degree.
- **G4 — The registry's own `QuadratureSpec` docstring contradicts its own
  entries.** `registry.py:259-261` says "*Lebedev and level-symmetric rules are
  **not** [half_range_clean] (they are designed to span the full sphere)*", while
  `registry.py:469` sets `half_range_clean=True` for `LevelSymmetricSN`.
  `[M]` **The entries are right, the docstring is wrong**:
  ```
    lebedev(5)             w(z>0)/w_tot = 0.366666666666667   equator nodes: 4
    lebedev(7)             w(z>0)/w_tot = 0.328571428571429   equator nodes: 8
    lebedev(11)            w(z>0)/w_tot = 0.429453262786596   equator nodes: 8
    level_symmetric(4)     w(z>0)/w_tot = 0.500000000000000   equator nodes: 0
    level_symmetric(8)     w(z>0)/w_tot = 0.500000000000000   equator nodes: 0
    product(4,8)           w(z>0)/w_tot = 0.500000000000000   equator nodes: 0
  ```
- **G5 — `axis_aligned` means different things per family.** `GaussLegendre1D`
  is tagged `axis_aligned=True` with the comment "*1-D nodes ARE the axis
  (μ-axis)*" — a category confusion: its nodes live on an interval, not on `S²`,
  so the flag is not comparable with Lebedev's (`[M]` 6 genuinely on-axis nodes).
  The `ProductQuadrature` comment "*depends on phi grid; conservatively False*" is
  half right — `[M]` it depends on the **parity of `n_mu`** (odd `n_mu` gives GL a
  node at `μ=0`) *and* the `φ` grid — and the exact-vs-tolerance split (1 exact vs
  4 at `1e-14`) is #325's defect showing through the flag.
- **G6 — No test measures any rule's actual exactness above degree 2** for the
  sphere families. `tests/numerics/test_rules_sphere.py` stops at
  `test_*_integrates_constants` (deg 0) and `test_*_quadratic_*` (deg 2), i.e.
  inside the sector `O_h` invariance makes free. That is why F1 survived.
- **G7 — No invariant is tested ACROSS realizations.** The one place the project
  has ratified cross-realization testing is #123's ≥2-quadrature protocol, and it
  is scoped to CP rank-N closures only. There is no
  "*every rule in the registry satisfies `∫ p dμ = ∫ p dx` for all `p` in its
  declared space*" parametrized test — which is exactly the test that ≥2
  realizations would let you write, and exactly the test that would have caught F1.
- **G8 — Not verified, flagged.** Whether any production/test configuration pairs
  `Quadrature.level_symmetric` with `scattering_order ≥ 2` at a single call site.
  Grep shows file-level co-occurrence in ~20 test modules; I did not resolve it to
  call sites. See §7.3.
- **G9 — Out of scope but noted.** `mc/solver.py:410` samples
  `theta = np.pi * rng.random()` — uniform in `θ`, not in `cos θ`. Not filed per
  the brief.

---

## 10. Reproduction

Probes (read-only, in the job tmp dir, not committed):
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe_exactness.py`,
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe_selector.py`.
Run with `.venv/bin/python`. All `[M]` blocks above are verbatim stdout.

---

## 11. Close re-check (L-007 / L-012)

Re-ran `git status --short` as the last action. One file appeared **during** this
dispatch that was not present at open:

- `.claude/plans/quadrature_machinery_campaign.md` *(untracked, created mid-dispatch)*

I did not read it — this survey is an independent inventory taken against the
tree, not a reconciliation against a plan. If the campaign plan already covers any
finding here, treat the overlap as corroboration from a separate route. No
production file changed under me (`git diff --name-only` shows only `.claude/`
memory, the `vv-principles` skill, and `tests/sn/_test_helpers.py`, all present at
open).
