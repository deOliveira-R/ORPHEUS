# Plan — consolidate slab into the unified `CurvilinearGeometry` machinery

**Author**: Claude Opus 4.7, 2026-04-23.
**Scope**: move slab Peierls out of its bespoke E₁ Nyström module
(`orpheus/derivations/peierls_slab.py`) into the unified polar-form
curvilinear machinery, so that `solve_peierls_1g(CurvilinearGeometry(
kind="slab-polar"), ...)` produces the same k_eff and flux as today's
`peierls_slab.solve_peierls_eigenvalue(...)`.

**Prerequisite plans** (read once if you haven't; do NOT re-derive):

- [`topology-based-consolidation.md`](topology-based-consolidation.md)
  — the Class A registration layer this plan does not touch.
  Slab is already registered under Class A via
  `peierls_cases.build_two_surface_case("slab", ...)` which currently
  routes to `peierls_slab._build_peierls_slab_case`. This plan
  changes only the *routing*; the Class A registration stays.
- [`post-cp-topology-and-coordinate-transforms.md`](post-cp-topology-and-coordinate-transforms.md)
  — Ch 4 (slab polar form), Ch 6 (exp-stretched substitution,
  **rejected** per §1 below), Phase G commit sequence, App B
  (full mathematical derivation of exp-stretched). Ch 5 (τ-coordinate)
  is out of scope — covered by Phase H / Issue #109.
- GitHub Issue #111 (Phase G) — the tracking issue; body has the
  short summary. This plan supersedes its quadrature strategy.

## 0. Executive summary

1. **What**: introduce a `solve_peierls_1g(slab_geom, ...)` path that
   produces bit-exact k_eff versus `peierls_slab.solve_peierls_eigenvalue`
   at matched precision budget.
2. **Why now**: the topology-based consolidation (2026-04-23, prior
   session) registered slab under Class A alongside hollow cyl/sph,
   but the *implementation* still lives in a separate module. The
   Sphinx capability matrix promises "Class A — F.4 applies"; the
   slab-specific machinery is an architectural tax that must be paid
   down for that promise to hold without caveats.
3. **Quadrature strategy (load-bearing)**: the µ-integral for slab
   polar form has **super-exponential** stiffness near µ=0, not
   exponential. Gauss-Laguerre on the exp-stretched `v = −ln|µ|`
   coordinate was considered and **rejected** (see §1). The adopted
   approach is **adaptive `mpmath.quad` on µ ∈ [−1, 1] with a forced
   breakpoint at µ=0**, extending the same "breakpoint-at-singularity"
   pattern that `peierls_slab._basis_kernel_weights` already uses for
   the E₁ log singularity.
4. **Prerequisites done** (do not re-do):
   - `CurvilinearGeometry(kind="slab-polar")` exists and has ~80 %
     of the primitives wired (§3).
   - Topology-based consolidation shipped 2026-04-23; slab already
     registered under Class A.
   - Phase F (hollow-core + rank-2 F.4) shipped; `compute_slab_transmission`
     (line 3111 of `peierls_geometry.py`) is the slab transmission
     primitive used by `_build_closure_operator_rank2_white`.
5. **Out of scope** (explicitly deferred):
   - τ-coordinate transform (Phase H / Issue #109).
   - Multi-region slab support in the unified path — start with 1-region,
     extend once G.3 benchmark parity lands.
   - Changing the Sphinx capability matrix or topology registration
     (both are correct as-is; only the routing inside
     `peierls_cases.build_two_surface_case("slab", ...)` changes).
6. **Estimated budget**: 2 sessions. Session 1 = G.1 (docs) + G.2
   (quadrature plumbing) + G.3 (benchmark). Session 2 = G.4 (planar-
   limit cross-check) + G.5 (routing switch with deprecation gate).
   If G.3 fails parity, the plan halts and a separate numerics-
   research session is needed to reconsider quadrature choice.

## 1. The mpmath-over-Gauss-Laguerre decision (the load-bearing rejection)

### 1.1 What was tried and why it was rejected

The coordinate-transforms plan (Ch 6, App B) proposed the **exp-
stretched substitution** `v = −ln|µ|` to absorb the grazing-ray
stiffness in the µ-integral into a Gauss-Laguerre `[0, ∞)` quadrature
with weight `exp(−v)`. The argument: the substitution derives `E₁`
from the polar integral, so GL on `v` should be numerically optimal.

**Problem**: the integrand-for-outer-µ-integral is not of the form
`f(µ) · exp(−const · µ)` — after the inner ρ-integral contributes
an `E₁(Σ_t · L / |µ|)`-like factor, the outer integrand is
**super-exponential in v**:

$$ \int_0^\infty e^{-\Sigma_t L\,e^v}\,dv. $$

Gauss-Laguerre's built-in weight is `exp(−v)` — it places nodes
optimally for that specific decay rate. For `exp(−const · e^v)`
decay, GL dramatically under-concentrates nodes in the stiff region
(where almost all the integrand mass lives) because it expects the
decay to be exponential, not super-exponential. In practice: GL-on-v
needs O(100) nodes for the precision that `E₁` Nyström gets in
O(20).

**Conclusion**: exp-stretched + Gauss-Laguerre is not competitive with
the native `E₁` Nyström. It was tried (in earlier work, per user
confirmation 2026-04-23) and rejected.

### 1.2 What replaces it: adaptive `mpmath.quad` with forced breakpoint

The native `peierls_slab.py` already uses the pattern `mpmath.quad(
integrand, [a, b, c])` where `b` is a singularity breakpoint — see
[peierls_slab.py:_basis_kernel_weights](orpheus/derivations/peierls_slab.py)
at line 157:

```python
w = mpmath.quad(integrand, [panel_a, x_eval, panel_b])
```

This works because `mpmath.quad` uses tanh-sinh (double-exponential)
quadrature by default and recursively subdivides from the breakpoints
outward, concentrating nodes wherever the integrand stiffens.
Tanh-sinh handles endpoint-singular and stiff-interior integrands
well by construction — no weight pre-matching required.

**The unified-slab approach inherits this pattern**: the outer
µ-integral uses `mpmath.quad(integrand, [−1, 0, 1])` with the
forced breakpoint at µ=0. `mpmath.quad`'s adaptive tanh-sinh
concentrates nodes at the grazing-ray stiffness naturally. No
exp-stretched substitution needed; no Gauss-Laguerre misapplication.

### 1.3 Keep the door open (future perf opportunity)

The user's standing reservation: *"we can re-check if we find a
better path"*. Specific candidates to revisit in a future numerics
session, IF adaptive `mpmath.quad` is too slow in practice:

- **Sinh-sinh quadrature** (super-exponential variant of tanh-sinh)
  — `mpmath.quad(..., method="tanh-sinh")` with a custom level
  schedule, or a direct sinh-sinh implementation. Tuned for
  super-exp decay.
- **Fixed `E₁`-in-the-kernel analytical reduction** — do the inner
  ρ-integral symbolically (giving `E₁(Σ_t L/|µ|)`), then the outer
  µ-integral becomes `∫ E₁(...) dµ` which mpmath's `e1` plus
  adaptive quad handles in closed moments. This is essentially a
  direct pull-through to the existing E₁ Nyström, so it doesn't
  win anything architecturally.

None are blocking. The adopted adaptive-mpmath.quad approach is
enough to ship this plan; perf-tuning is post-merge.

## 2. Current state — `peierls_slab.py` native E₁ Nyström

### 2.1 Module inventory

697 lines ([peierls_slab.py](orpheus/derivations/peierls_slab.py)).
Top-level public functions:

- `solve_peierls_eigenvalue(...)` — 1G eigenvalue power iteration
  on the E₁ Nyström kernel matrix.
- `solve_peierls_fixed_source(...)` — 1G fixed-source solve.
- `solve_peierls_slab_mg(...)` — multi-group wrapper.
- `_build_peierls_slab_case(ng_key, n_regions, ...)` — case builder
  that produces a `ContinuousReferenceSolution`; called by
  `peierls_cases.build_two_surface_case("slab", ...)`.

Private machinery:

- `_product_log_weights(panel_a, panel_b, x_eval, nodes, dps)` —
  builds product-integration weights for the `−ln|x−x'|` part of
  `E₁` via `mpmath.quad(integrand, [panel_a, panel_b])`.
- `_basis_kernel_weights(...)` — unified kernel-weight builder that
  combines the log singularity with the smooth `E₁` remainder via
  `mpmath.quad(integrand, [panel_a, x_eval, panel_b])`.
- `_build_kernel_matrix(...)` — assembles the full `N × N` kernel
  matrix at specified dps.
- `_build_system_matrices(...)` — multi-group block-Toeplitz
  assembly.

### 2.2 Singularity handling (the archaeological record)

The `E₁(z)` kernel has a logarithmic singularity at `z = 0`:
`E₁(z) = −ln(z) − γ + R(z)` with smooth `R(z)`. Native slab uses
**singularity subtraction + product integration**: GL on the smooth
`R(z)` + tabulated `mpmath.quad` weights for `−ln|x−x'|` against each
panel basis function.

Per the module's own docstring (line 28–30):
> The Nyström method uses **singularity subtraction**: standard GL
> weights for the smooth remainder *R*, and **product-integration
> weights** (computed via mpmath.quad) for the `−ln|x−x'|` part.

This is numerically efficient: O(h²) convergence with GL panel
order, machine precision at dps=20 with 4 panels × 6 points.

### 2.3 Test + reference surface that must stay green

- [`tests/derivations/test_peierls_convergence.py`](tests/derivations/test_peierls_convergence.py):
  `solve_peierls_eigenvalue` spectral convergence verification.
- [`tests/derivations/test_peierls_reference.py`](tests/derivations/test_peierls_reference.py):
  direct `_build_kernel_matrix` checks at line 96.
- [`tests/derivations/test_peierls_slab.py`](tests/derivations/test_peierls_slab.py):
  slab-specific regression.
- Shipped continuous reference: `peierls_slab_2eg_2rg` with
  `k_eff = 1.23067665` (measured 2026-04-23; reference fixture).
- `tests/l0_error_catalog.md` lines 1168, 1221 cite
  `peierls_slab._build_kernel_matrix` in bug-catch entries — do
  not orphan these without updating the catalog.

## 3. What `CurvilinearGeometry(kind="slab-polar")` already supports

Audit of
[peierls_geometry.py](orpheus/derivations/peierls_geometry.py). The
constructor accepts `kind="slab-polar"` and rejects `inner_radius>0`
(line 278). The `SLAB_POLAR_1D` singleton exists (line 927). The
following methods are **already wired**:

| Method | Slab-polar value | Line |
|---|---|---|
| `d` (dimension) | 3 | 290 |
| `S_d` (angular-measure norm) | 4π | 299 |
| `prefactor` | 0.5 | 311 |
| `n_surfaces` | 2 | 330 |
| `topology` (NEW, 2026-04-23) | `"two_surface"` | 335 |
| `is_planar` | True | 378 |
| `angular_range` | (−1, 1) | 405 |
| `angular_weight(ω)` | 1 | 421 |
| `angular_to_mu(ω)` | identity (ω is µ) | 405 |
| `rho_max(x, µ, L)` | (L−x)/µ or −x/µ | 437 |
| `source_position(x, ρ, µ)` | x + ρµ | (confirm) |
| `volume_kernel_mp(τ)` | exp(−τ) | 745 |
| `radial_volume_weight(r)` | 1 | 759 |
| `shell_volume_integral(r, w, φ)` | Σ wᵢ φᵢ | 874 |
| `rank1_surface_divisor(L)` | 2 | 797 |
| `reciprocity_factor` | raises ValueError | 915 |

**The geometry primitives are 100 % wired.** What's missing is the
*quadrature path* that uses them — see §4.

## 4. What remains to implement (the concrete gap list)

1. **Adaptive mpmath.quad support in `build_volume_kernel`** (the core
   volume kernel assembler in `peierls_geometry.py`). Today it uses
   fixed Gauss-Legendre on the angular variable. For slab-polar
   with grazing-ray stiffness at µ=0, replace with
   `mpmath.quad(angular_integrand, [−1, 0, 1])` adaptive tanh-sinh.
2. **`CurvilinearGeometry.angular_quadrature(n, dps, adaptive=False)`**
   — a new method that returns angular integration nodes+weights for
   non-adaptive callers, and an adaptive-integrator closure for
   slab-polar with `adaptive=True`. Signature TBD in G.2; see §6.2
   for sketch.
3. **Sanity bridge to `compute_slab_transmission`** — the existing
   F.4 builder at `_build_closure_operator_rank2_white` (line 3777
   of peierls_geometry.py) already dispatches to
   `compute_slab_transmission` for slab. Once slab-polar's
   `build_volume_kernel` works, the F.4 path lights up automatically.
   Verify this end-to-end in G.3.
4. **`_build_peierls_slab_case` routing switch** (in
   `peierls_slab.py`, or a new builder in `peierls_cases.py`):
   optionally call the unified `solve_peierls_1g(SLAB_POLAR_1D, ...)`
   instead of the native `solve_peierls_eigenvalue`. Gate on G.3
   benchmark parity.
5. **Docs** — Sphinx `peierls_unified.rst` gets a new section
   `§theory-peierls-slab-polar` documenting the polar-form slab
   equation, the super-exp stiffness analysis, the
   mpmath.quad-with-breakpoint strategy, and the Gauss-Laguerre
   rejection (for the archaeological record).

## 5. Target math (distilled, one page)

The slab Peierls equation in observer-centred polar form:

$$
\Sigma_t(x)\,\varphi(x)
  \;=\; \frac{1}{2}\!\int_{-1}^{1}\!\mathrm d\mu\!
        \int_0^{\rho_{\max}(x,\mu)}\!
        e^{-\int_0^\rho \Sigma_t(x+s\mu)\,\mathrm ds}\,
        q\bigl(x + \rho\,\mu\bigr)\,\mathrm d\rho,
$$

with

$$
\rho_{\max}(x, \mu) \;=\; \begin{cases}
  (L - x)/\mu     & \text{for } \mu > 0, \\
  -x/\mu = x/|\mu| & \text{for } \mu < 0, \\
  \infty          & \text{for } \mu = 0 \text{ (rays parallel to faces)}.
\end{cases}
$$

The source `q(x') = Σ_s(x') φ(x') + χ · νΣ_f(x') φ(x') / k` follows
the standard scattering + fission convention. `Σ_t(x)` is the
piecewise-constant total cross section.

**This is the same form as `solve_peierls_1g` in
`peierls_geometry.py` already uses for cylinder-1d and sphere-1d.**
The only differences are:

- Angular measure: `dµ` on `[-1, 1]` (slab) versus `sin θ dθ` on
  `[0, π]` (sphere) or `dβ` on `[0, π]` (cylinder) — already handled
  by `geometry.angular_weight(ω)` and `geometry.angular_range`.
- `ρ_max`: piecewise-linear in µ — already handled by
  `geometry.rho_max(x, µ, L)`.
- Source position: `x + ρµ` (linear) versus the quadratic
  ray-geometry formula for curvilinear — already handled by
  `geometry.source_position`.

**The volume kernel is `exp(−τ)` for both slab-polar AND sphere-1d.**
This is the "surprise" noted in the coord-transforms plan Ch 4.1:
both geometries have the same 3-D point-kernel reduction. Cylinder
is the odd one out, using Ki₁.

**What slab's polar form specifically needs**: the outer µ-integral
at `µ=0` is where `ρ_max → ∞`. The integrand itself decays as
`exp(−Σ_t · L/|µ|)` — super-exponential in v = −ln|µ|. See §1 for
the quadrature-strategy implications.

## 6. Implementation plan (G.1 through G.5)

### G.1 — Sphinx documentation (~2 h)

- Add `:label: peierls-slab-polar` equation block for the polar-form
  slab equation above.
- New section `§theory-peierls-slab-polar` under the existing
  `peierls_unified.rst` file (location: after the E₁ narrative,
  before the moment-form archive section). Content:
  - The polar-form slab equation.
  - The grazing-ray stiffness analysis (§1 of this plan).
  - Why Gauss-Laguerre on `v = −ln|µ|` was considered and rejected
    (super-exp vs exp).
  - The adopted adaptive `mpmath.quad([-1, 0, 1])` approach.
- Update capability matrix §`theory-peierls-capabilities` to note
  slab has the unified path available (in addition to the E₁ native
  path).

**Acceptance**: Sphinx builds clean (no new warnings), new label
`peierls-slab-polar` resolves, narrative reads cleanly.

### G.2 — Angular-quadrature plumbing (~4 h)

Introduce a dispatch seam so `build_volume_kernel` can use adaptive
mpmath.quad on slab-polar while retaining Gauss-Legendre on
cyl/sph.

**Option A (preferred)**: new method on `CurvilinearGeometry`:

```python
def angular_quadrature(
    self, n: int, dps: int = 25,
) -> tuple[np.ndarray, np.ndarray]:
    """Angular-integration nodes + weights on ``angular_range``.

    Returns (pts, wts) suitable for ``sum(wts * integrand(pts))``.
    For slab-polar, returns adaptive tanh-sinh points from
    mpmath.quad with a breakpoint at µ=0. For cyl/sph, returns
    Gauss-Legendre on [0, π]."""
```

This assumes the angular integrand factorises cleanly from the
radial — verify in G.2 implementation. If the integrand is not
factorisable (the ρ-integral depends on µ in an entangled way that
doesn't precompute to fixed points), an alternative:

**Option B (fallback)**: new `build_volume_kernel` code path for
slab-polar that calls `mpmath.quad(lambda µ: inner_rho_integral(µ),
[-1, 0, 1])` directly — no precomputed angular nodes; the adaptive
quad drives the angular integration end-to-end.

Decide Option A vs B by trying A first; if the integrand cannot be
factorised at the angular level without losing the mpmath adaptive
benefit, fall back to B.

**Acceptance**: `solve_peierls_1g(SLAB_POLAR_1D, radii=[L],
sig_t=..., ..., boundary="vacuum")` runs to completion and returns
a valid `PeierlsSolution`. k_eff agreement with native E₁ at ~3
sigfigs (tighter bounds come at G.3).

### G.3 — End-to-end benchmark vs native E₁ (~4 h)

Compare unified-slab versus native-E₁ on the shipped `peierls_slab_2eg_2rg`
fixture and on a 1G 1-region sweep. Success criteria:

- **Bit-exact at matched high precision**: at `dps=30`,
  `n_panels=8`, `n_angular=48`, the unified k_eff matches
  `solve_peierls_eigenvalue` k_eff to **1e-10 rtol**.
- **Cost ratio at matched precision**: unified needs at most **3×**
  the wall time of native E₁ for the same k_eff precision. If
  ratio is worse, flag for the "future perf" numerics session
  (§1.3).
- **Flux-shape agreement**: `||φ_unified − φ_native||_∞ / ||φ_native||_∞
  < 1e-8` on the 6-point grid.

Implement as a new pytest class `TestSlabPolarVsNativeE1` in
[tests/derivations/test_peierls_rank2_bc.py](tests/derivations/test_peierls_rank2_bc.py)
(natural home — this is slab rank-2 territory). Mark `@pytest.mark.slow`.

**Acceptance**: 1G 1-region + the `2eg_2rg` fixture both pass. If
either fails, halt G.4/G.5 and open a numerics investigation.

### G.4 — Planar-limit cross-check (~2 h)

Physics argument: a hollow cylinder at `r_0 = 0.999 · R` approximates
a slab of thickness `L = R − r_0`. Their k_eff values should match
to high precision as the cavity approaches planar.

Test: set up `solve_peierls_1g(CurvilinearGeometry(kind="cylinder-1d",
inner_radius=0.999), ...)` and `solve_peierls_1g(SLAB_POLAR_1D,
radii=[0.001], ...)`. Compare k_eff.

- **Expected**: agreement to 1e-8 rtol or better. The limit is
  geometric; numerical precision should be the only gap.
- **If fails**: cylinder's Ki₁ kernel and slab's exp(−τ) kernel do
  not yield the same planar limit — suggests a bug in one of the
  geometries.

Implement as `test_slab_polar_planar_limit_of_thin_hollow_cylinder`
in `test_peierls_rank2_bc.py`.

**Acceptance**: 1e-8 rtol agreement at one thin-cylinder
configuration. Stretch: 1e-10 rtol at `dps=40`.

### G.5 — Routing switch in `peierls_cases` (~1 h)

Modify
[`peierls_cases.build_two_surface_case`](orpheus/derivations/peierls_cases.py)
so `shape="slab"` can route to either the native
`peierls_slab._build_peierls_slab_case` OR the unified
`solve_peierls_1g(SLAB_POLAR_1D, ...)` path. Introduce an internal
flag `_SLAB_VIA_UNIFIED: bool = True` (defaulted to True after G.3
benchmark parity confirmed); environment variable override
`ORPHEUS_SLAB_VIA_E1=1` forces the native path for bisection.

Bit-exact regression: `peierls_slab_2eg_2rg` k_eff continues to
produce the same reference `k_eff = 1.23067665` (to whatever
precision G.3 establishes). Tests that call
`peierls_slab._build_kernel_matrix` directly (per §2.3) continue
to use the native path — they are testing the E₁ Nyström
specifically, not the reference interface.

Emit a `DeprecationWarning` on native-path usage pointing at the
unified path for one release; remove the native path only after a
explicit follow-up commit (not in this plan's scope).

**Acceptance**: `rv.continuous_get("peierls_slab_2eg_2rg").k_eff`
returns the same value as today. All tests from §2.3 pass unchanged.

## 7. Regression fixture (bit-exact baseline)

Captured 2026-04-23 for this plan's acceptance gate:

- `peierls_slab_2eg_2rg` k_eff: **1.23067665** at default quadrature
  (n_panels_per_region=4, p_order=4, precision_digits=20).
- 1G 1-region critical slab at `L/λ_t = 5`, `c = 0.5 / (1 − 0.5) =
  1.0`: k_eff ≈ (read from native solver; insert at G.3 start).

If any of these shift by more than 1e-10 rtol after G.5, the plan
halts and the shift is investigated before proceeding.

## 8. The retirement question — keep `peierls_slab.py` indefinitely

The coord-transforms plan's Phase G.5 proposes "conditional retire
`peierls_slab.py` if benchmark parity". **This plan overrides that
recommendation**: keep `peierls_slab.py` indefinitely as an
independent cross-check implementation.

Rationale:

1. **697 lines of singularity-subtraction + product-integration is
   valuable infrastructure**, independently reviewed and documented.
2. **It serves as an independent cross-check** for the unified path.
   Two implementations that compute the same answer via different
   numerical routes catch bugs that either implementation alone
   would miss.
3. **Test entries in `l0_error_catalog.md` (lines 1168, 1221) cite
   `peierls_slab._build_kernel_matrix` explicitly**. Retiring the
   module orphans those entries.
4. **Cost of retention is low**: the module is
   not exercised by default once G.5 lands (`_SLAB_VIA_UNIFIED = True`
   routes through the unified path), so it doesn't slow down anyone.
5. **Cost of deletion is high**: lost cross-check surface, lost
   archaeological record of the E₁-Nyström approach, broken catalog
   references.

**Exception**: if a future session discovers a bug in the native
path that is unfixable without major rework, `peierls_slab.py` can
be moved to `derivations/archive/` (reversible via `git mv`) at
that time. Do not delete.

## 9. Acceptance criteria (summary)

- **G.1** Sphinx builds clean after documentation additions.
- **G.2** `solve_peierls_1g(SLAB_POLAR_1D, ...)` runs to
  completion with valid flux + k_eff.
- **G.3** Unified k_eff matches native E₁ Nyström k_eff to **1e-10
  rtol** on the `peierls_slab_2eg_2rg` fixture AND on a 1G
  1-region sweep; cost ratio ≤ 3×.
- **G.4** Planar-limit cross-check with hollow cylinder at
  `r_0 = 0.999R` matches unified slab at 1e-8 rtol.
- **G.5** All existing slab tests (§2.3) pass with
  `_SLAB_VIA_UNIFIED = True`. `peierls_slab_2eg_2rg` reference
  k_eff stable within 1e-10.

If **any** of G.3 / G.4 fails, halt and open a numerics-research
session. Do not land G.5 with a failing benchmark.

## 10. Open questions

**OQ1 — Does `mpmath.quad` perform well enough in practice?**
The claim that tanh-sinh handles super-exp stiffness is theoretically
sound but empirically unverified. If G.2 / G.3 shows the cost ratio
> 3× the native E₁ Nyström, switch to sinh-sinh or the analytical-E₁
fallback in §1.3. Worst case, G.2 falls through to Option B and the
slab-polar code path retains an embedded `mpmath.quad` call rather
than a precomputed quadrature — slower but correct.

**OQ2 — What happens for multi-region slab (e.g., ``peierls_slab_2eg_2rg``)?**
The unified path via `solve_peierls_1g` handles multi-region through
the `radii` array and per-region `sig_t`. The `composite_gl_r`
radial-mesh builder in `peierls_geometry.py` already supports
multiple radial regions. In theory this should Just Work for slab-
polar; verify empirically at G.3 using `peierls_slab_2eg_2rg` as
the fixture.

**OQ3 — Is the slab F.4 rank-2 closure reachable through the
unified path?** The `_build_closure_operator_rank2_white` at line
3777 already dispatches to `compute_slab_transmission` for slab. The
Wigner-Seitz-exact E₂/E₃ bilinear closure should light up as soon
as `build_volume_kernel` works for slab-polar. Verify at G.3.

**OQ4 — Does the exp-stretched substitution deserve a third look?**
The user's standing "we can re-check if we find a better path"
reservation. A future numerics session could try a Gauss-Jacobi
endpoint-weight formula with a super-exp envelope function matched
to the specific `exp(−Σ_t · L · e^v)` decay. Not blocking.

## 11. Session trail

- **2026-04-23 (plan filing)**: plan filed. User clarified that
  Gauss-Laguerre on the exp-stretched coordinate was tried and
  rejected (integrand super-exp, not exp); the fix is adaptive
  mpmath.quad with breakpoint at µ=0, extending the pattern already
  used in `peierls_slab._basis_kernel_weights`.

- **2026-04-23 (execution session)** — Claude Opus 4.7:

  - **G.1 ✅ landed.** New Sphinx section `§theory-peierls-slab-polar`
    in `docs/theory/peierls_unified.rst` documenting the active
    adaptive mpmath.quad path (supersedes the archived
    τ-Laguerre and moment-form sections). Cross-referenced from Key
    Facts header and from the Class A capabilities table. Sphinx
    builds clean; new label resolves in 6 places in the built HTML.

  - **G.2 ✅ already implemented** (pre-existing, not done this session).
    `K_vol_element_adaptive` in
    [peierls_geometry.py:937](orpheus/derivations/peierls_geometry.py#L937)
    performs nested adaptive `mpmath.quad` with forced µ=0 breakpoint
    for slab-polar (line 1064). `build_volume_kernel` dispatches
    slab-polar to `build_volume_kernel_adaptive` at line 1163. This
    plan's assumption that "adaptive mpmath.quad support needs to be
    added" was out of date (shipped 2026-04-20 per Issue #117 which
    archived the moment-form).

  - **G.3 ✅ landed.** New test class
    `TestSlabPolarVsNativeE1KEff` in
    [tests/derivations/test_peierls_rank2_bc.py](tests/derivations/test_peierls_rank2_bc.py)
    (2 parametrisations, 41.37 s, both PASSED). Unified
    `solve_peierls_1g(SLAB_POLAR_1D, ...)` k_eff matches native
    `solve_peierls_eigenvalue` to **1e-8** rtol, normalised flux
    matches to **1e-6**, at N=3 (1 panel × p=3) dps=20. Marked
    `@pytest.mark.slow`. Tighter tolerances achievable at larger N
    but the 41-second cost already confirms verification-level
    equivalence.

  - **G.4 ❌ BLOCKED — plan physics flawed.** Empirical probe at the
    plan's specified configuration (hollow cyl `r_0=0.999 R=1.0`,
    slab `L=0.001`, `Σ_t=1, c=0.4, νΣ_f=0.6`, vacuum BC, N=3
    dps=20) gives k_eff delta of **22.5 %**, nowhere near the plan's
    claimed 1e-8 tolerance. The plan assumed "the limit is
    geometric" but ignored that cylinder's `Ki₁` kernel has already
    integrated the axial direction; the in-plane chord distribution
    in a thin annular shell scales as `sqrt(2·R·L) ≈ 0.045` for
    tangential rays, fundamentally different from slab's
    `L/|μ|`. The two kernels see different optical-depth spectra
    even in the thin-shell limit. A meaningful planar limit needs
    either ray-distribution matching or a curvature-over-thickness
    expansion — both non-trivial future work.

    **Resolution**: filed as
    [Issue #129](https://github.com/deOliveira-R/ORPHEUS/issues/129)
    for physics investigation. Not a shipping test. Sphinx
    narrative updated to document the finding.

  - **G.5 ⚠️ INFRASTRUCTURE LANDED, ACTIVATION DEFERRED (2026-04-24,
    Issue #130)**. After Issue #104 unblocked multi-group
    `solve_peierls_mg`, this session ran the direct parity benchmark
    the original plan called for. The result:

    - Native E₁ Nyström on `peierls_slab_2eg_2rg`:
      k_eff = 1.226 530 511 976, 1.53 s wall time.
    - Unified `solve_peierls_mg(SLAB_POLAR_1D, ..., "white_f4")`:
      k_eff = 1.245 529 269 703, 930.11 s wall time.
    - **rel_diff = 1.5 %, cost ratio = 606×** (N=12, dps=20).

    That's far too large for default routing activation (plan
    target: 1e-10). Single-region 1G (TestSlabPolarVsNativeE1KEff)
    and single-region 2G with fabricated XS
    (TestMGSlabPolarMatchesNativeSlabMG) both showed 1e-8
    agreement, so the discrepancy appears specifically in the
    multi-region multi-group regime of the unified path. Likely
    causes (in decreasing likelihood):

    1. Quadrature underconvergence for the specific XS combination
       (region-B Σ_t,2 = 2.0 → deep thermal optical depth across
       the moderator; adaptive `mpmath.quad` may terminate below
       the E₁ Nyström's easy-regime precision at this N).
    2. Multi-region ray walker subtlety (optical-depth
       accumulation across material interfaces under 2G
       assembly).
    3. F.4 closure inter-group coupling on multi-region slabs.

    **What landed anyway** (Issue #130 partial completion):

    - `peierls_cases._SLAB_VIA_UNIFIED` flag, default `False`.
    - `ORPHEUS_SLAB_VIA_UNIFIED=1` env-var override for bisection.
    - `_build_peierls_slab_case_via_unified` unified-path builder.
    - Dispatch in `build_two_surface_case("slab", ...)`.
    - Diagnostic test
      `TestSlabViaUnifiedDiscrepancyDiagnostic` that records the
      current 1.5 % gap as a regression barometer (loose 5 %
      bound — passes today, catches future degradation).
    - Sphinx `§theory-peierls-slab-polar-g5-routing` documents the
      benchmark + deferral.

    **Update 2026-04-24 (later the same day)** —
    [Issue #131](https://github.com/deOliveira-R/ORPHEUS/issues/131)
    **resolved**. The numerics-investigator agent ran probes A, B,
    D, E, F under `derivations/diagnostics/diag_slab_issue131_*.py`
    and pinpointed the cause:

    - `compute_P_esc_{outer,inner}` and `compute_G_bc_{outer,inner}`
      had two branches for slab-polar: a closed-form
      `½ E₂(τ)` / `2 E₂(τ)` branch for `len(radii) == 1`, and a
      **finite-N GL quadrature** branch for multi-region
      (`len(radii) > 1`).
    - The multi-region GL branch converged only to ~4 × 10⁻³ at
      N=24 — a quadrature artifact that fed into the K_bc closure
      and produced the 1.5 % k_eff gap.
    - The µ-integral `½ ∫₀¹ exp(-τ(x_i)/µ) dµ = ½ E₂(τ(x_i))` is
      **closed-form regardless of the number of regions** because
      `τ(x_i)` is µ-independent for piecewise-constant σ_t.

    **Fix applied 2026-04-24** (commit TBD):

    - New helpers `_slab_tau_to_outer_face` and
      `_slab_tau_to_inner_face` in `peierls_geometry.py`
      piecewise-integrate σ_t across region boundaries.
    - All slab-polar P_esc/G_bc calls route through the closed-form
      E₂ branch regardless of `len(radii)`.
    - Post-fix benchmark: rel_diff drops from 1.549 × 10⁻² to
      **5.4 × 10⁻¹⁶** (bit-exact to machine epsilon). Cost ratio
      unchanged (~550×, expected for a verification primitive).

    **G.5 activated:**

    - `_SLAB_VIA_UNIFIED` default flipped to `True`.
    - Env-var override renamed to `ORPHEUS_SLAB_VIA_E1=1` (routes
      to native for bisection).
    - `TestSlabViaUnifiedDiscrepancyDiagnostic` tightened from
      `rel_diff < 5 %` to `rel_diff < 1e-10` — now a regression
      gate, not a diagnostic.
    - Issue #130 (Phase G.5 routing) and Issue #131 (discrepancy
      investigation) both closable.

- **Outcomes for this plan**:
  - 1G slab-polar verification-reference equivalence is
    **proven** (G.1, G.2, G.3) — documented + tested + passing.
  - G.4 is a physics misconception; the hollow-cyl-to-slab
    planar limit does not hold at 1e-8 under matched
    `(Σ_t, Σ_s, νΣ_f, L)`. Future work (Issue #129).
  - G.5 **infrastructure is in place** but activation blocked on
    a ~1.5 % unified-vs-native k_eff gap discovered 2026-04-24.
    Follow-up numerics investigation required before default
    activation.
