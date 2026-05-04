# Post-CP Plan — Topological unification & coordinate-transformation exploits

**Status**: DEFERRED — to be executed AFTER the CP module unification campaign currently running in parallel (`.claude/plans/next-session-phase-4-3-cp-unification.md`, Phases A–E). The current in-flight work covers:

- Phase A: sphere Peierls (Phase 4.3)
- Phase B: flat-source CP unification (cp_geometry.py, BickleyTables retirement)
- Phase C: solver.py cleanup (#102)
- Phase D: higher-rank white BC, multi-group, V&V audit, perf
- Phase E: wrap-up

**This plan picks up after E.** It captures mathematical insights identified in session-N that would be lost without explicit documentation. The insights fall into two families:

1. **Topological unification** (Part I): observation that 2-boundary geometries — slab, hollow cylinder, hollow sphere — share a common topological class, distinct from solid curvilinear (1-boundary + coordinate singularity). Extending ORPHEUS to hollow-core geometries and then absorbing slab into the unified curvilinear framework as the planar limit.

2. **Coordinate-transformation exploits** (Part II): several classical coordinate transforms that move numerical difficulty into places where our quadrature handles it gracefully. The `τ`-as-integration-variable transform (Monte Carlo delta-tracking for Nyström) is the most underexploited and stands to pay the biggest dividends.

**Audience**: future Claude Code sessions picking up after the CP campaign completes. The plan is written with the full mathematical context that session-N has in its head; future sessions will enter cold and this document is their primer.

---

## 0. Context primer — what session-N knows that you don't

Before executing any phase in this plan, read:

1. `docs/theory/peierls_nystrom.rst` (Sessions N and N+1 wrote this — ~1500+ lines)
2. `orpheus/derivations/peierls_geometry.py` — the `CurvilinearGeometry` abstraction
3. `orpheus/derivations/peierls_cylinder.py` — the thin facade over `peierls_geometry`
4. `orpheus/derivations/peierls_sphere.py` — added in Phase A (if not yet done, this plan is premature)
5. The companion plan `next-session-phase-4-3-cp-unification.md` (in `~/.claude/plans/` — that plan's container is ephemeral; the salient excerpts are copied to Appendix A of THIS document for redundancy)

### 0.1 The unified Peierls polar form

At the end of the CP campaign, every curvilinear 1-D radial geometry is an instance of:

$$\Sigma_t(r)\,\varphi(r) \;=\; \frac{\Sigma_t(r)}{S_d}\!
  \int_{\Omega_d}\!d\Omega\!\int_0^{\rho_{\max}(r,\Omega)}\!
    \kappa_d(\Sigma_t\rho)\,q(r'(\rho,\Omega,r))\,d\rho
  + S_{\rm bc}(r)$$

with `CurvilinearGeometry` holding `kind ∈ {"cylinder-1d", "sphere-1d"}` and supplying the four geometry-specific primitives: `angular_weight(Ω)`, `rho_max(r, Ω, R)`, `source_position(r, ρ, Ω)`, `volume_kernel(τ)` (i.e. `Ki_1` or `e^{-τ}`).

### 0.2 The 3-tier kernel integration hierarchy

The kernels E_n / Ki_n / exp are successive integrations of a single 3-D point kernel. The hierarchy (documented in `peierls_nystrom.rst`):

```
Level 0 (point):             e^{-τ}/(4πR²)
     ↓ integrate over symmetry (2/1/0 dims)
Level 1 (pointwise Peierls):  E_1/2   Ki_1/(2π)   e^{-τ}/(4π)
     ↓ integrate over source region (flat-source average)
Level 2 (partial current):    E_2     Ki_2          e^{-τ}
     ↓ integrate over target region (flat-source average)
Level 3 (flat-source CP):     E_3     Ki_3          e^{-τ}
```

With `E_n' = -E_{n-1}`, `Ki_n' = -Ki_{n-1}`. This identity drives much of what follows: the τ-coordinate transform in Part II exploits the FACT that the first-level integrand e^{-τ} is geometry-invariant; higher levels inherit this uniformity.

### 0.3 What's NOT unified yet at post-CP entry

- **Slab** remains in its own module (`peierls_slab.py`) using the E_1 Nyström with singularity subtraction. Not unified with curvilinear because E_1 has a log singularity at τ=0 that the standard curvilinear polar form doesn't have.
- **Hollow-core geometries** (annular cylinder with r_0 > 0 inner radius, spherical shell with r_0 > 0 inner radius) are NOT supported. The `radii` array in current code is outer-radii of concentric annuli starting from r=0.
- **τ-coordinate transform** is not used anywhere. Ray integration is in ρ (geometric distance), which forces the kernel to track Σ_t(r(ρ)) through multi-region walkers.

This plan addresses all three.

---

## Part I — Topological unification

### Chapter 1: The 2-boundary topological class

#### 1.1 Classifying 1-D radial geometries

Every 1-D radial geometry in ORPHEUS falls into one of two classes based on topology:

| Class | Number of boundaries | Number of coordinate singularities | Members |
|-------|---------------------|-----------------------------------|---------|
| **A** (2-boundary) | 2 (inner + outer) | 0 | slab `[0, L]`; hollow cylinder `[r_0, R]` with `r_0 > 0`; hollow sphere `[r_0, R]` with `r_0 > 0` |
| **B** (1-boundary + singularity) | 1 (outer) | 1 (at `r=0`) | solid cylinder `[0, R]`; solid sphere `[0, R]` |

**Key observation**: all Class-A geometries are diffeomorphic to the interval `[a, b]`. The slab `[0, L]` and the hollow cylinder `[r_0, R]` and the hollow sphere `[r_0, R]` are topologically the same space — just with different metrics.

Class-B geometries have an interior coordinate singularity at `r=0` that is **geometrically regular** — rays pass through `r=0` without incident — but represents a coordinate artifact, not a physical boundary.

#### 1.2 BC structural consequences

The 2-boundary topology has immediate consequences for boundary-condition closures:

- **Class A** geometries have TWO boundaries, each carrying its own scalar `J^+` and `J^-` partial currents. White BC gives rank-2 closure (the slab's `E_2`-based outer product that `peierls_slab.py` already implements). Reflective BC also rank-2.
- **Class B** geometries have ONE boundary (the outer surface). White BC for radially-symmetric flux collapses to rank-1 at the flat-source level (our current Phase 4.2 implementation). At the pointwise Nyström level, rank-1 is the Mark/isotropic approximation and carries 21% error at R=1 MFP (Issue #100 + Phase 4.2 C7 cylinder limitation).

**Adding hollow-core support lifts Class B into Class A**. The hollow cylinder / sphere inherits the slab's rank-2 BC structure. This is structurally clean.

#### 1.3 Why this distinction matters for us

The refactor in Phase B (of the parallel plan) built `cp_geometry.py` as a unified flat-source CP infrastructure. Class B was unified first (cyl + sph in one framework). Slab remains apart because Class-A-specific BC machinery differs.

**If ORPHEUS adds hollow-core support**, every 1-D radial geometry becomes Class A, and the slab's BC machinery generalizes over the full set. This is a very clean unification.

### Chapter 2: Hollow-core geometry implementation

#### 2.1 API changes to `CurvilinearGeometry`

Current `CurvilinearGeometry` (post-CP-campaign):

```python
@dataclass(frozen=True)
class CurvilinearGeometry:
    kind: str  # "cylinder-1d" or "sphere-1d"
    # implicit: inner boundary at r=0 (solid)
```

Extension:

```python
@dataclass(frozen=True)
class CurvilinearGeometry:
    kind: str  # "cylinder-1d", "sphere-1d", "cylinder-hollow", "sphere-hollow"
    inner_radius: float = 0.0  # r_0. Defaults to 0 (solid case).
    # When inner_radius > 0, the geometry is hollow.
```

#### 2.2 Ray geometry changes

For hollow cylinder/sphere, a ray from observer `r_i` in direction `Ω` may hit either the OUTER boundary at `r=R` or the INNER boundary at `r=r_0`, whichever comes first.

**Outer intersection** (as in current code):

$$\rho_{\max}^{\rm outer}(r_i, \Omega) = -r_i \cos\Omega + \sqrt{r_i^2 \cos^2\Omega + R^2 - r_i^2}$$

**Inner intersection** (new): solve `r_i² + 2 r_i ρ cosΩ + ρ² = r_0²`:

Discriminant: `Δ = r_i² cos²Ω + r_0² - r_i²`. For `Δ < 0`, the ray misses the inner cavity entirely. For `Δ ≥ 0`, the two roots are:

$$\rho^{\rm inner}_\pm = -r_i \cos\Omega \pm \sqrt{r_i^2 \cos^2\Omega + r_0^2 - r_i^2}$$

Both roots are non-negative when the ray enters and exits the cavity. The ray enters at `ρ^{\rm inner}_-` (closer) and exits at `ρ^{\rm inner}_+` (farther).

**Physical interpretation** (depends on the inner BC):

- **Vacuum inner BC**: the ray terminates at `ρ^{\rm inner}_-`. `ρ_{\max}(r_i, \Omega) = min(ρ^{\rm outer}, ρ^{\rm inner}_-)`.
- **Reflective inner BC**: the ray bounces off at `ρ^{\rm inner}_-`; for Nyström purposes this is equivalent to the "mirror" continuation.
- **Transmissive (hollow cavity, no material)**: the ray traverses the cavity with zero Σ_t, then continues through the cavity and re-enters the material on the other side at `ρ^{\rm inner}_+`. Continues until it hits the outer surface at some `ρ > ρ^{\rm inner}_+`.

The **transmissive case is the most important** for fuel-pin applications (central plenum / coolant channel). It requires the optical-depth walker to handle the cavity as a zero-Σ_t segment.

#### 2.3 Optical-depth walker extension

Current `optical_depth_along_ray` walks annular crossings where `Σ_t` is piecewise constant on `[0, r_1], [r_1, r_2], ..., [r_{N-1}, R]`. For hollow with cavity, the walker treats the cavity interval `[−ρ^{\rm inner}_+, −ρ^{\rm inner}_-]` (or its mirror) as zero-Σ_t. This is a small patch to the existing walker.

The key insight: when the τ-coordinate transform (Part II, Chapter 5) is in place, the cavity becomes a ρ-jump discontinuity in the `ρ(τ)` map but contributes zero to τ itself. In τ-space, the cavity is **invisible** — it's just a segment of the ρ-axis where `dρ/dτ = ∞` (Σ_t = 0) but `τ` doesn't advance.

#### 2.4 Inner BC rank-1 closure

For a hollow geometry with inner radius `r_0`, the inner surface has area:
- Cylinder: `A_{\rm inner} = 2π r_0` (per unit z)
- Sphere: `A_{\rm inner} = 4π r_0²`

White BC at the inner surface mirrors the outer rank-1 closure but with separate `J^+_{\rm in}`, `J^-_{\rm in}` scalars:

$$J^-_{\rm in} = J^+_{\rm in}, \qquad \text{and separately} \qquad J^-_{\rm out} = J^+_{\rm out}$$

Giving **rank-2** total closure. The structure mirrors slab's `E_2` outer product with two face contributions.

Need new geometry-specific helpers:
- `compute_P_esc_inner(r_i, radii, sig_t)` — probability of escape to the inner surface
- `compute_G_bc_inner(r_i, ...)` — surface-to-volume Green's function from inner surface

Both analogous to the outer versions but integrated over the inner surface geometry.

#### 2.5 Test configurations

Real physics test cases for hollow-core:

- **BWR-style fuel pin** with central plenum: `r_0 = 0.1`, `R = 0.5`, annular UO₂ between (simplified: single-region annulus of fissile material).
- **TRISO fuel kernel**: `r_0 = 0`, `R = 0.1`, ..., but with multiple shells (coated particle). Can represent as hollow annular layers.
- **Pebble-bed coolant channel**: sphere with central cavity.

Verification: compare against Sanchez-type tabulations for annular cylinders (if findable via Zotero/lit-researcher after the server comes back).

### Chapter 3: The planar limit connecting slab to hollow curvilinear

#### 3.1 Mathematical setup

Consider a hollow cylinder with `r_0 → ∞` and `R = r_0 + L` fixed. Physical intuition: as the inner radius grows, the annular region becomes a "thin strip" with local curvature `κ = 1/r_0 → 0`. In the limit of zero curvature, this is a slab of thickness L.

#### 3.2 Asymptotic expansion of ray geometry

For observer at `r_i = r_0 + x` where `x ∈ [0, L]`, in direction `Ω = β`:

$$r'(\rho, \beta, r_0 + x) = \sqrt{(r_0+x)^2 + 2(r_0+x)\rho\cos\beta + \rho^2}$$

$$= r_0\sqrt{\left(1 + \frac{x}{r_0}\right)^2 + \frac{2(1+x/r_0)\rho\cos\beta}{r_0} + \frac{\rho^2}{r_0^2}}$$

For `r_0 \gg x, \rho`, Taylor expand:

$$r' = r_0 + x + \rho\cos\beta + O(1/r_0)$$

So `r' - r_0 = x + ρ cosβ + O(1/r_0)`. The position in the "slab coordinate" `x' = r' - r_0` is just `x + ρ cosβ` — the slab's linear source formula. 

Similarly, `ρ_max`:

$$\rho_{\max}^{\rm outer}(r_i, \beta) = -r_i\cos\beta + \sqrt{r_i^2\cos^2\beta + R^2 - r_i^2}$$

For `r_i = r_0 + x` and `R = r_0 + L`:

$$R^2 - r_i^2 = (r_0+L)^2 - (r_0+x)^2 = 2r_0(L-x) + (L^2 - x^2)$$

To leading order in `1/r_0`:

$$\rho_{\max}^{\rm outer} \to -(r_0+x)\cos\beta + \sqrt{(r_0+x)^2\cos^2\beta + 2 r_0(L-x) + O(1)}$$

For `cos β > 0` (outward rays): `ρ_{\max} → (L-x)/cosβ`.

For `cos β < 0` (inward rays toward inner surface `r=r_0`), the outer-boundary intersection goes off to larger `ρ` values; the inner-boundary intersection gives `ρ^{\rm inner}_- → x/|cosβ|`.

These are **exactly the slab ρ_max formulas** with `μ = cos β`. The planar limit reduces hollow cylinder to slab.

#### 3.3 Kernel asymptotics: Ki_1 → E_1 / 2 in the limit

The Ki_1 kernel of the cylinder came from z-integration of the 3-D point kernel over the (formally) infinite z-axis. In the planar limit, the z-axis is still infinite but now the geometric structure is "slab-like". The Ki_1 → E_1/2 reduction happens when we REDO the angular integration.

Explicitly: for the slab, integrating the 3-D point kernel over transverse directions (the 2-D plane perpendicular to x), then over the remaining angular direction μ, gives `E_1(τ)/2`. For the cylinder, the procedure is z-first (giving Ki_1) then 2-D polar (giving the polar form). In the planar limit, these two procedures must give the same answer (both are just dim-3 → dim-1 in different orders).

The identity:

$$\int_0^\infty \mathrm{Ki}_1(\tau\sqrt{1+u^2})\,du = \frac{1}{2}E_1(\tau) + \text{curvature corrections}$$

where `u` is a transverse-plane coordinate in the planar-limit expansion. (This identity would be valuable to verify numerically during Phase G.2.)

#### 3.4 Verification: use planar limit as a regression check

During Phase G (slab absorbed into unified framework), we can verify that the hollow-cylinder solution at `r_0 = 100` (say, 100 MFP) reproduces the slab solution of equal thickness to high precision. This is a non-trivial cross-check that exercises both the hollow-core code (Chapter 2) AND the planar-limit correctness.

### Chapter 4: Slab as an instance of unified framework

#### 4.1 Direct approach (preferred)

Don't go through hollow curvilinear; just use slab's own observer-centered polar form with `μ = cos(angle from x-axis)`:

$$\varphi(x) = \frac{1}{2}\!\int_{-1}^1\!d\mu\!\int_0^{\rho_{\max}(x, \mu)}\!e^{-\Sigma_t \rho}\,q(x + \rho\mu)\,d\rho$$

with `ρ_max(x, μ>0) = (L-x)/μ` and `ρ_max(x, μ<0) = x/|μ|`.

This fits the `CurvilinearGeometry` API with:
- `kind = "slab"`
- `angular_weight = 1` (uniform in μ)
- `angular_range = (-1, 1)`
- `rho_max(x, μ)` = the piecewise formula above
- `source_position(x, ρ, μ) = x + ρ μ`
- `volume_kernel(τ) = e^{-τ}` (same as sphere!)
- `prefactor = 1/2`

**Slab and sphere share the `e^{-τ}` kernel in the polar form.** (Surprise — but obvious in retrospect; both are 3-D problems fully characterized by the 3-D point kernel.)

#### 4.2 The grazing-ray stiffness

The slab polar form has `ρ_max = L/|μ| → ∞` as `μ → 0`. The integral converges (exp(-Σ_t L/|μ|) → 0 super-exponentially) but the INTEGRAND is stiff near `μ = 0`: most of the contribution comes from a vanishingly small range of μ.

**Quadrature options**:

1. **Standard Gauss-Legendre on μ ∈ [-1, 1]**: fails — `O(1/N²)` convergence at best due to stiffness.
2. **Gauss-Jacobi with endpoint weight**: `(1-μ²)^α` weight could help but doesn't absorb the specific exp-stiff behavior.
3. **Exp-stretched coordinate** (Chapter 6): `v = -ln|μ|`, so `μ = e^{-v}` with `v ∈ [0, ∞)`. Then Gauss-Laguerre on `v` handles the exp behavior naturally.
4. **Double exponential (Tanh-Sinh)**: very robust for endpoint-singular integrands; standard in mpmath.

Option 3 is the mathematically elegant choice — it's literally the substitution that derives E_1 from the polar integral. Option 4 is the practical default if option 3 is cumbersome to implement.

#### 4.3 Retire peierls_slab.py or keep as specialist form?

Two options:

**Option α — Fully retire peierls_slab.py**. The unified `CurvilinearGeometry(kind="slab")` subsumes it. Pros: one module for all geometries; cleaner architecture. Cons: the native E_1 Nyström in peierls_slab.py uses singularity subtraction + product integration which is numerically efficient; the unified polar form may need MORE quadrature points to achieve the same precision (stiff grazing rays).

**Option β — Keep peierls_slab.py as specialist, make unified framework support slab too**. Pros: no loss of the efficient E_1 implementation; gives users a choice. Cons: two ways to do the same problem; users must choose; code duplication.

**Recommendation: Option α**, subject to benchmarking. If the unified slab is within 2x the E_1 version's precision-per-quadrature-point, retire the E_1 version. If it's slower, keep both with documented trade-offs.

---

## Part II — Coordinate-transformation exploits

### Chapter 5: The τ-coordinate transform (most underexploited)

#### 5.1 Motivation

In Nyström ray integration, we compute:

$$\int_0^{\rho_{\max}} e^{-\tau(\rho)} q(r'(\rho)) d\rho$$

where `τ(ρ) = ∫_0^ρ Σ_t(r(s)) ds`. For multi-region, `Σ_t(r(s))` is piecewise constant with jumps at annular crossings.

The current code evaluates `e^{-τ(ρ)}` explicitly — requiring the walker to compute the piecewise integral `τ(ρ)` at each quadrature node.

**The τ-coordinate transform** uses τ as the integration variable:

$$\tau \equiv \int_0^\rho \Sigma_t(r(s))\,ds, \qquad d\tau = \Sigma_t(r(\rho))\,d\rho$$

Change of variable:

$$\int_0^{\rho_{\max}} e^{-\tau(\rho)}\,q(r'(\rho))\,d\rho \;=\; \int_0^{\tau_{\max}}\!e^{-\tau}\,\frac{q(r'(\tau))}{\Sigma_t(r'(\tau))}\,d\tau$$

where `τ_max = τ(ρ_max)` and `r'(τ) = r'(ρ(τ))` with `ρ(τ)` the inverse map.

#### 5.2 Why this is so nice

1. **The kernel `e^{-τ}` is geometry-INDEPENDENT and medium-INDEPENDENT.** It's always just `exp(-τ)` on `τ ∈ [0, τ_max]`. Multi-region structure, vacuum cavities, anisotropic scattering — all absorbed into the `ρ(τ)` map and the `1/Σ_t` Jacobian.

2. **Quadrature efficiency.** Gauss-Laguerre on `[0, ∞)` with `e^{-τ}` weight IS THE IDEAL quadrature for this integrand. Every node is placed optimally. Compare: currently we use Gauss-Legendre on `[0, ρ_max]`, which over-samples regions where `e^{-τ}` is already tiny.

3. **Hollow cavities are TRIVIAL.** In a ρ-integral, traversing a cavity contributes `∫_{ρ_a}^{ρ_b} e^{-τ(ρ_a)} q(r'(ρ)) dρ` (τ doesn't advance in the cavity since Σ_t = 0). In τ-space, the cavity is just a SINGLE POINT — `τ = τ(ρ_a) = τ(ρ_b)` — and the integration skips over it automatically.

4. **Multi-region uniformity.** No piecewise-constant Σ_t in the kernel. The map `ρ(τ)` is piecewise linear with different slopes in each region; the kernel never sees this.

#### 5.3 Monte Carlo connection

The τ-coordinate is **literally** the sampling variable used in Monte Carlo delta-tracking / Woodcock tracking. In MC, particles sample their next collision distance as `τ = -ln(ξ)` where `ξ ∈ (0, 1]` is uniform, then walk the ρ-coordinate until they reach that τ. The inverse map `ρ(τ)` is computed ON THE FLY by the ray walker.

For Nyström, the same idea: compute the ρ(τ) map once (per ray), then use Gauss-Laguerre in τ.

MC literature that's relevant:
- Woodcock et al. 1965 (original delta-tracking)
- Brown & Martin 2003 (delta-tracking review)
- Leppänen 2010 (Serpent code implementation details)

#### 5.4 Implementation plan

Extend `CurvilinearGeometry` with a τ-as-ρ option:

```python
def optical_depth_along_ray_with_map(self, r_obs, cos_omega, tau_max,
                                      radii, sig_t):
    """Return list of (τ_k, ρ_k) breakpoints for the piecewise-linear ρ(τ) map.
    Used by τ-coordinate quadrature to compute ρ(τ) at arbitrary τ."""
    # ... walker that returns the piecewise-linear map instead of just final τ
```

Add a `use_tau_coordinate: bool` flag to `build_volume_kernel`:

```python
def build_volume_kernel(geometry, r_nodes, ..., use_tau_coordinate=False):
    if use_tau_coordinate:
        # Use Gauss-Laguerre in τ
        tau_pts, tau_wts = gauss_laguerre(n_rho, dps)  # ∫_0^∞ e^-τ f(τ) dτ
        # τ_max bound: truncate Laguerre at τ_max via endpoint-handling
        # Map back to ρ via piecewise-linear ρ(τ)
        ...
    else:
        # Existing GL-in-ρ code
        ...
```

#### 5.5 Expected benefits

- **Optically thick cells**: 2-4x fewer quadrature points for same precision (exp-weighted Gauss is optimal)
- **Optically thin cells**: ≈ same (both methods converge quickly)
- **Multi-region**: CLEANER code, no piecewise-in-kernel logic
- **Hollow cavities** (when Phase F lands): natural, no special-case code for cavity traversal

#### 5.6 Caveats and risks

- The map `ρ(τ)` requires tracking piecewise-linear breakpoints at annular boundaries. For `Σ_t = 0` cavities, `dρ/dτ = ∞` creating a singularity in the map. Handle by making the cavity a "τ-jump" (no contribution) rather than smooth.
- Gauss-Laguerre on `[0, ∞)` but our domain is `[0, τ_max]`. Truncate via complementary error function or use tanh-sinh on the finite interval.
- Re-benchmarking required after implementing — the classical Gauss-Legendre might still win for optically thin problems or when α `n_rho` is already large.

#### 5.7 Phase-H commit sequence

1. **H.1**: Derive the τ-transform mathematically and add to `peierls_nystrom.rst` as a new section "Coordinate transformations for Nyström quadrature".
2. **H.2**: Implement `optical_depth_along_ray_with_map` returning the piecewise-linear ρ(τ) breakpoints.
3. **H.3**: Add `use_tau_coordinate=True` code path to `build_volume_kernel`. Gate on benchmark.
4. **H.4**: Benchmark: precision vs quadrature-node-count for both ρ- and τ-integration on canonical test cases (thin cylinder, thick cylinder, multi-region, hollow).
5. **H.5**: If τ wins, make it the default. Else keep it optional for optically-thick cases.

### Chapter 6: Exp-stretched coordinate for slab

#### 6.1 The slab grazing-ray problem

In observer-centered polar form for slab:

$$\varphi(x) = \frac{1}{2}\!\int_{-1}^1\!d\mu\!\int_0^{\rho_{\max}(x,\mu)}\!e^{-\Sigma_t\rho}\,q(x+\rho\mu)\,d\rho$$

with `ρ_max = (L-x)/μ` (for μ>0) or `x/|μ|` (for μ<0). As `μ → 0`, `ρ_max → ∞`. The integrand has a stiff exp-decay structure.

#### 6.2 The exp-stretched substitution

Define `v = -ln|μ|` for `μ ∈ (0, 1]`, so `μ = e^{-v}`, `dμ = -e^{-v} dv`, `1/μ = e^v`:

$$\int_0^1 e^{-\Sigma_t L/\mu}/\mu\,d\mu \;=\; \int_0^\infty e^{-\Sigma_t L e^v}\,dv$$

(Check for slab-specific context: if the inner integrand is not just `e^{-Σ_t L/μ}` but involves `q(x + ρμ)` with the ρ-integral done first — then the substitution acts on the `1/μ`-stiff factor.)

The integrand `e^{-Σ_t L e^v}` decays super-exponentially for large `v`, making Gauss-Laguerre in `v` ideal. Specifically, if we substitute once more `u = Σ_t L e^v`:

$$\int_0^\infty e^{-\Sigma_t L e^v}\,dv = \int_{\Sigma_t L}^\infty e^{-u} \frac{du}{u} = E_1(\Sigma_t L)$$

This recovers the E_1 form of the slab problem. The exp-stretched `v` is the coordinate that implicitly generates `E_1` — so the polar-form slab with this substitution is MATHEMATICALLY EQUIVALENT to the E_1 Nyström, just evaluated differently.

#### 6.3 Why this is useful

It lets us put slab in the unified `CurvilinearGeometry` framework without giving up the E_1-equivalent numerical efficiency. The slab quadrature becomes:

1. Integrate over `v` using Gauss-Laguerre (absorbs `e^{-v}` weight naturally).
2. For each `v`, compute `μ = e^{-v}` and the ρ-integral on `[0, L/|μ|]`.
3. The exp-stiff behavior at `μ → 0` becomes uniform behavior at `v → ∞`, handled naturally by GL.

Or, if we do the ρ-integral analytically first (giving E_1-dependence), then the outer integral is over `v` with a well-behaved integrand.

#### 6.4 Implementation sketch

```python
class CurvilinearGeometry:
    ...
    def slab_angular_quadrature(self, n, dps=25):
        """Gauss-Laguerre on v ∈ [0, ∞) with μ = e^(-v) substitution.
        Returns (v_k, w_k, μ_k, w_μ_k) for both directions.
        """
        if self.kind != "slab":
            raise ValueError("Only for slab")
        # Gauss-Laguerre on [0, ∞) with weight e^-v
        v_pts, v_wts = mpmath.gauss_quadrature(n, "laguerre")
        mu_pts = np.array([float(mpmath.exp(-v)) for v in v_pts])
        # The weight transformation: integral over μ = integral over v with factor μ = e^-v
        # which is already the Laguerre weight, so w_μ = w_v (up to normalization)
        mu_wts = np.array([float(w) for w in v_wts])
        # Return both positive-μ and negative-μ branches
        return mu_pts, mu_wts, -mu_pts, mu_wts
```

Then `build_volume_kernel` for slab uses this quadrature instead of standard GL.

### Chapter 7: Gauss-Jacobi endpoint weights (fallback method)

#### 7.1 The general method

For integrals `∫_a^b (x-a)^α (b-x)^β f(x)\,dx` with `α, β > -1`, Gauss-Jacobi quadrature has weight function `(x-a)^α (b-x)^β` built in. For polynomial `f` of degree `≤ 2n-1`, the rule is exact. For smooth `f`, convergence is exponential (for bounded domain).

#### 7.2 Application to chord-form `1/√(r'²-y²)`

The chord-form integrand has `1/√(r'²-y²) = 1/√((r'-y)(r'+y))`. On `r' ∈ [y, R]`, the singular endpoint is `r' = y` with weight `(r'-y)^{-1/2}`. Gauss-Jacobi with `α = -1/2, β = 0` absorbs this.

**But we already solved this by pivoting to polar form in Phase 4.2.** Gauss-Jacobi is a fallback if we ever need the chord form again (e.g., for connecting to closed-form CP formulas).

#### 7.3 Application to E_1 log singularity

The slab E_1 kernel has `E_1(x) ~ -ln(x) - γ` near `x=0`. No simple Gauss-Jacobi weight captures the log behavior directly; the existing singularity-subtraction approach in `peierls_slab.py` is the right method. Noted for completeness.

#### 7.4 Recommendation

Document Gauss-Jacobi availability in the code but don't default to it. Use when a specific problem has a known power-law singularity at an endpoint.

### Chapter 8: Legendre angular expansion for higher-rank BC

(This chapter extends the parallel-plan's Phase D.1 with additional theoretical detail. Read D.1 first; this chapter adds mathematical context that session-N had in mind but didn't write down.)

#### 8.1 The rank-1 failure mode (recap)

Current Mark/isotropic white-BC closure sets `J^-(Ω) = J^-/π` (isotropic inward over the hemisphere). In reality, the angular distribution of re-entering particles IS anisotropic — the transport solution tells us ψ(Ω) is concentrated in directions close to the inward normal (because those paths to the surface are shorter and less attenuated).

**Rank-1** = truncating the Legendre expansion of ψ(Ω) at `n=0` (isotropic). Higher orders (n ≥ 1) capture the anisotropy.

#### 8.2 Angular expansion basis

The natural basis for white-BC in curvilinear geometry is **Legendre polynomials in `μ_s = Ω · n̂_s`** where `n̂_s` is the inward surface normal:

$$\psi(\Omega) = \sum_{n=0}^{N} (2n+1) \,c_n\, P_n(\mu_s)$$

Orthogonality on `μ_s ∈ [0, 1]` (inward hemisphere):

$$\int_0^1 P_n(\mu_s) P_m(\mu_s) d\mu_s = \frac{\delta_{nm}}{2n+1} \,\times\, \text{(convention-dependent factor)}$$

#### 8.3 Rank-N closure

Each `c_n` is a separate scalar unknown. The outgoing angular flux has its own Legendre expansion:

$$\psi_{\rm out}(\Omega) = \sum_n (2n+1) d_n P_n(\mu_{s,{\rm out}})$$

White BC per mode: `d_n = c_n` (each mode reflects to itself).

The closure is rank-N: N_modes scalar unknowns per surface (x2 for 2-boundary geometries).

#### 8.4 Computing `d_n` from interior source

For each mode `n`, `d_n = (2n+1)/ ∫_{\rm outward} P_n(\mu_s) · ψ_{\rm out}(\Omega)\,dΩ`, where `ψ_{\rm out}` is constructed by attenuated transport from the interior source:

$$\psi_{\rm out}(\Omega) = \int_0^{\rho_{\max}(x_s, \Omega_{\rm back})} e^{-\tau(\rho)} q(r'(\rho, \Omega_{\rm back})) d\rho$$

where `Ω_back` is the reversed direction (from surface inward toward source) and `x_s` is the specific surface point.

Inner product against `P_n(μ_s)` gives `d_n` — one scalar per mode per surface.

#### 8.5 Convergence rate

The Legendre expansion converges **geometrically** for smooth `ψ(Ω)` and **polynomially** for `ψ` with a boundary-layer structure. The physical boundary layer at a surface under white BC is typically O(1 MFP) wide, so the angular flux has a discontinuity at the grazing angle (μ_s = 0).

Empirical guidance:
- **Thick cells (R > 5 MFP)**: N = 1 (Mark, current code) gives < 1% error.
- **Moderate cells (R ~ 2-5 MFP)**: N = 3-5 should give < 0.1% error.
- **Thin cells (R < 2 MFP)**: N = 8-16 may be needed.

**Benchmark plan**: implement for N ∈ {1, 2, 4, 8} and measure convergence vs Sanchez tie-point at multiple R values.

#### 8.6 Cost

Rank-N closure adds N × N_surfaces scalar unknowns per group. For 1D radial single-surface: N scalar unknowns. For the linear-algebra cost, we solve a rank-N perturbation of the volume kernel — the Sherman-Morrison formula extends to rank-N (Woodbury identity) to avoid rebuilding the full matrix. Cost per eigenvalue iteration: O(N·N_r) + O(N³) vs the rank-1 case's O(N_r) + O(1).

For N ≤ 16 and N_r ~ 50, this is trivially fast.

### Chapter 9: Davison's u = r·φ (historical; keeps for connection to analytic results)

#### 9.1 The classical sphere trick

For spherically-symmetric flux on `[0, R]`, define `u(r) = r · φ(r)`. The Peierls equation for sphere transforms into a 1-D integral equation on u:

$$\Sigma_t(r) u(r) = \int_0^R K(r, r') u(r') dr' + S_{\rm bc}(r)$$

with `u(0) = 0` as a natural BC (since `φ(0)` is finite, `u(0) = 0·φ(0) = 0` exactly).

#### 9.2 Why session-N didn't use it

The observer-centered polar form NEVER needs to evaluate `φ(0)` directly. Rays from an interior observer at `r_i > 0` pass through `r=0` uneventfully (it's just a regular interior point geometrically). The Lagrange basis on composite panels doesn't include `r=0` as a node (panels start with GL nodes interior to `[0, r_1]`). So the Davison substitution's benefit (regularizing `r=0`) doesn't materialize for polar-form Nyström.

#### 9.3 Where it IS useful

- **Analytical closed-form solutions**: Case-Zweifel's classical analytic solutions for bare spheres use the Davison substitution. If we ever want to compare against these at high precision, the substitution provides the numerical connection.
- **Differential-form transport**: if we ever want to write the Peierls equation in differential form (which we don't, currently), the Davison trick is needed to regularize the operator at r=0.
- **Asymptotic analysis**: the u = rφ form is the natural variable for the classical Milne problem and bare-sphere critical-radius asymptotic expansions.

#### 9.4 Recommendation

Don't implement. Document in `peierls_nystrom.rst` as a historical note. If ever needed, add as a specialist transform.

---

## Part III — Phase-by-phase implementation plan

### Phase F: Hollow-core geometry support

**Estimated effort**: 5-7 commits, ~1-2 days of focused work.

**Prerequisites**: 
- CP module unification campaign complete (Phases A–E from parallel plan)
- `CurvilinearGeometry` + `peierls_geometry.py` stable
- Sphere Peierls exists

**Commit sequence**:

- **F.1** `feat(derivations): add inner_radius to CurvilinearGeometry; update rho_max`
  - Extend dataclass with `inner_radius: float = 0.0`
  - Add inner-intersection root calculation in `rho_max`
  - Update `optical_depth_along_ray` to handle cavity segment
  - All existing tests (solid geometries) pass: `inner_radius = 0` is the default.

- **F.2** `feat(derivations): inner-surface P_esc_inner and G_bc_inner`
  - Implement `compute_P_esc_inner(geometry, ...)`
  - Implement `compute_G_bc_inner(geometry, ...)` for cylinder (sphere needs its own derivation, may defer)
  - Test: for solid geometry (`inner_radius = 0`), both return zeros (no inner surface).

- **F.3** `feat(derivations): rank-2 white-BC closure for hollow cylinder`
  - Extend `build_white_bc_correction` to handle two surfaces
  - Rank-2 structure: `K_bc = u_out⊗v_out + u_in⊗v_in`
  - Test: for solid geometry, reduces to rank-1 (current case).

- **F.4** `feat(derivations): hollow-sphere analogs of P_esc_inner, G_bc_inner`
  - Mirror F.2/F.3 for sphere
  - Add Sphinx derivation in `peierls_nystrom.rst` Chapter N "Hollow-core geometries"

- **F.5** `test(derivations): hollow-cylinder and hollow-sphere verification`
  - Annular cylinder with central cavity: critical radius vs published tables (if Zotero gives us data)
  - Spherical shell: analytical limits (thick shell → effectively infinite medium)

- **F.6** `docs(theory): hollow-core geometry section in peierls_nystrom.rst`
  - Full derivation of inner-surface machinery
  - Planar-limit connection to slab (briefly)
  - Physical motivation (fuel pins, pebble beds)

- **F.7** (optional) `feat(cp): hollow-core support in cp_geometry flat-source`
  - Mirror F.1-F.6 for flat-source CP
  - Needed if ORPHEUS wants to support hollow-fuel CP problems

### Phase G: Slab absorbed into unified framework

**Estimated effort**: 4-6 commits, ~1-2 days.

**Prerequisites**: Phase F (to validate via planar limit)

**Commit sequence**:

- **G.1** `docs(theory): coordinate transforms in Peierls Nyström (new section)`
  - New chapter in `peierls_nystrom.rst`: "Coordinate transformations for numerical stability"
  - Sections on τ-coordinate, exp-stretched μ, Gauss-Jacobi, Davison u=rφ
  - The content of THIS document's Part II, expanded

- **G.2** `feat(derivations): exp-stretched slab in CurvilinearGeometry`
  - Add `kind = "slab"` to `CurvilinearGeometry`
  - Implement `slab_angular_quadrature` using exp-stretched μ
  - Source position, ρ_max, volume_kernel specialized for slab

- **G.3** `test(derivations): slab Peierls in unified framework regression`
  - Verify unified-framework slab matches existing `peierls_slab` for all Phase-4.1 test cases
  - Precision target: ≤ 1e-6 relative difference
  - Speed target: ≤ 2x slower than E_1 Nyström at same precision

- **G.4** `test(derivations): planar-limit verification`
  - Solve hollow cylinder at r_0 = 100, R = 101, L = 1
  - Compare to slab at L = 1
  - Should agree to 1e-4 (curvature corrections ~ 1/r_0 ~ 1%)

- **G.5** (conditional) `refactor(derivations): retire peierls_slab.py`
  - Only if G.3 benchmark shows unified-framework slab is within 2x precision-per-node of E_1 form
  - Thin-facade peierls_slab.py delegating to `CurvilinearGeometry(kind="slab")`
  - All slab tests pass unchanged

- **G.6** `docs(theory): slab as unified framework instance (amendment)`
  - Update `peierls_nystrom.rst` to list slab among the unified geometries
  - Cross-reference the exp-stretched section

### Phase H: τ-coordinate retrofit

**Estimated effort**: 3-4 commits, ~1 day.

**Prerequisites**: None (can be done any time after the CP campaign). Orthogonal to Phases F, G.

**Commit sequence**:

- **H.1** `docs(theory): τ-coordinate as Nyström integration variable`
  - New section in `peierls_nystrom.rst` (or extend the coordinate-transforms chapter from G.1)
  - Full derivation, Monte Carlo connection, cost/benefit

- **H.2** `feat(derivations): optical_depth_along_ray_with_map`
  - Returns piecewise-linear ρ(τ) breakpoints for a ray
  - Testable in isolation

- **H.3** `feat(derivations): use_tau_coordinate=True option in build_volume_kernel`
  - Gauss-Laguerre on τ ∈ [0, τ_max] with endpoint truncation
  - Computes r'(τ) via breakpoints from H.2
  - Jacobian 1/Σ_t factored in

- **H.4** `test(derivations): τ-coordinate vs ρ-coordinate convergence comparison`
  - Same problem, both quadrature schemes, varying n_rho (or n_tau)
  - Measure: precision(k_eff) vs total wall time
  - Expected: τ wins for optically-thick cells (R ≥ 5 MFP); ρ wins for thin cells

- **H.5** (conditional) `feat: make τ-coordinate default for optically-thick cells`
  - Heuristic: if estimated τ_max > 5, use τ-coordinate
  - Else use ρ-coordinate
  - Dispatch hidden in `build_volume_kernel`

### Phase I: Cross-reference to higher-rank BC (main plan's Phase D.1)

Already planned in the parallel document. This plan's Chapter 8 adds theoretical detail that should be folded into Phase D.1 of that plan. Key contributions from Chapter 8:

- Legendre basis choice (expansion in `μ_s = Ω · n̂_s`, not azimuthal)
- Convergence-rate expectations for thin/moderate/thick cells
- Woodbury/Sherman-Morrison update for rank-N perturbations to avoid rebuilding the kernel matrix

When Phase D.1 is executed, the implementer should reference Chapter 8 of THIS document.

### Phase ordering summary

```
Post-CP campaign:
    (Phase D.1 — higher-rank BC, from main plan. Blocks some of this plan's phases.)

This plan:
    Phase F (hollow-core) ─┐
                           │
    Phase G (slab) ────────┤    ... independent of each other; H can be any-order
                           │
    Phase H (τ-coord) ─────┘
```

- F and G can be done independently, but F unlocks a nice regression check for G (planar limit).
- H is independent but orthogonal — could be done FIRST, and would benefit F's hollow-cavity computation (cavity becomes a τ-no-op).
- Optimal order: **H first** (provides τ-coordinate for all subsequent), then **F** (hollow-core using τ), then **G** (slab using polar form, benchmarking against F's planar limit).

---

## Part IV — Rich Sphinx documentation content

The user explicitly requested "excruciating detail" in the Sphinx docs for these topics. This section sketches the content for a dedicated Sphinx chapter to be added to `peierls_nystrom.rst` (or as a new file `docs/theory/peierls_nystrom_advanced.rst` if the primary page gets too long).

### 4.1 Proposed Sphinx structure

```
Part III — Beyond the basic unification (NEW chapter in peierls_nystrom.rst)

§20 Topology of 1-D radial geometries
   20.1 The 2-boundary class and the 1-boundary+singularity class
   20.2 Diffeomorphism between Class-A members
   20.3 BC structural consequences (rank-2 vs rank-1)
   20.4 Hollow-core as Class-A extension

§21 The planar limit: slab as hollow-cylinder → ∞
   21.1 Curvature expansion
   21.2 Ki_1(τ) → E_1(τ)/2 as r_0 → ∞
   21.3 Numerical verification strategy
   21.4 Why direct polar-form slab is preferred in practice

§22 Coordinate transformations in Nyström quadrature
   22.1 Principles: transforms don't eliminate singularities, they relocate them
   22.2 Jacobian-absorbing: s² = r'² - y² (already used)
   22.3 Optical-depth coordinate τ = ∫Σ_t ds
   22.4 Exp-stretched μ = e^(-v) for slab grazing rays
   22.5 Gauss-Jacobi endpoint weights
   22.6 Davison u = r·φ (historical)

§23 Higher-rank white-BC: Legendre angular expansion
   23.1 The Mark/isotropic truncation and its error
   23.2 Angular expansion basis (μ_s vs Ω)
   23.3 Rank-N closure structure (Woodbury update)
   23.4 Convergence rates by cell thickness

§24 Monte Carlo connections
   24.1 Delta-tracking ≡ τ-coordinate sampling
   24.2 Next-event estimator ≡ fixed-ray Nyström
   24.3 Where deterministic and stochastic methods diverge
```

### 4.2 Equation labels to introduce

All with `:vv-status: documented`:

- `peierls-tau-coordinate-transform` — the variable change τ(ρ)
- `peierls-planar-limit-cylinder-to-slab` — the r_0 → ∞ asymptotic
- `peierls-exp-stretched-mu` — the μ = e^{-v} substitution for slab
- `peierls-rank-n-bc-closure` — higher-rank white BC
- `peierls-hollow-core-rho-max` — ray with two boundary intersections
- `peierls-delta-tracking-equivalence` — MC delta-tracking ↔ τ-Nyström

### 4.3 Key numerical evidence to accumulate (as sections are added)

- Thick-vs-thin convergence comparison table: Gauss-Legendre-in-ρ vs Gauss-Laguerre-in-τ at various cell thicknesses
- Rank-N white-BC convergence table: k_eff error vs N at R ∈ {1, 2, 5, 10} MFP
- Planar-limit table: hollow-cylinder k_eff at varying r_0 → slab k_eff
- Slab polar-form vs E_1-Nyström: precision per quadrature node

### 4.4 Citations to add

- **Woodcock et al. 1965** — original delta-tracking paper (MC connection)
- **Brown & Martin 2003** — delta-tracking review (via OpenAlex; DOI searchable)
- **Leppänen 2010** (Serpent paper) — practical delta-tracking implementation
- **Bell & Glasstone 1970** — already cited; Section 2.7 has the Davison substitution
- **Case, de Hoffmann & Placzek 1953** — classical bare-sphere critical radii using Davison

---

## Part V — Mathematical appendix (the derivations that must not be lost)

### App A: τ-coordinate transformation — full derivation

Starting from the polar-form Peierls integrand for a ray from observer `r_i` at angle `Ω`:

$$I(r_i, \Omega) = \int_0^{\rho_{\max}(r_i, \Omega)} e^{-\tau(\rho)} q(r'(\rho)) d\rho$$

where `τ(ρ) = ∫_0^ρ Σ_t(r(s)) ds`.

**Step 1**: Change variable `τ = τ(ρ)`. Since `τ(ρ)` is monotonically non-decreasing (Σ_t ≥ 0), it has an inverse `ρ = ρ(τ)` (multi-valued at cavity Σ_t=0 segments, but we pick the canonical branch).

**Step 2**: Compute differentials. `dτ = Σ_t(r(ρ)) dρ`, so `dρ = dτ / Σ_t(r(ρ))`. At `ρ(τ)`:

$$dρ = \frac{1}{\Sigma_t(r(ρ(τ)))} dτ = \frac{1}{\Sigma_t(r'(τ))} dτ$$

where `r'(τ) = r(ρ(τ))` is the source position as a function of τ.

**Step 3**: Substitute:

$$I = \int_0^{\tau_{\max}} e^{-\tau}\,q(r'(\tau))\,\frac{d\tau}{\Sigma_t(r'(\tau))}$$

with `τ_max = τ(ρ_max)`.

**Step 4**: For cavity segments where `Σ_t = 0`, `τ` does not advance (`dτ = 0`) while `ρ` does. These segments contribute nothing to the integral (in the τ-coordinate they are single points). Numerically: the cavity interval `[ρ_a, ρ_b]` maps to `τ(ρ_a) = τ(ρ_b)`, a point of measure zero in τ-space.

**Step 5**: For multi-region with `Σ_t` piecewise constant, `ρ(τ)` is piecewise linear with slope `1/Σ_{t,k}` in each annular region. Representing `ρ(τ)` by its breakpoints `{(τ_k, ρ_k)}` at annular crossings allows exact evaluation at any τ.

**Quadrature**: Gauss-Laguerre on `[0, ∞)` with weight `e^{-τ}`, truncated at `τ_max`:

$$\int_0^{\tau_{\max}} e^{-\tau} g(\tau) d\tau \approx \sum_k w_k^{\rm Lag} g(\tau_k) \quad \text{for nodes } \tau_k < \tau_{\max}$$

Plus a correction for the tail `[τ_max, ∞)` which is subtracted (since `e^{-τ} g(τ)` is typically tiny there).

Alternative: Gauss-Legendre on `[0, τ_max]` directly. But Laguerre uses the exp weight natively.

### App B: Exp-stretched coordinate for slab — full derivation

The slab Peierls integrand (dropping the q(x+ρμ) factor temporarily to focus on the stiffness):

$$J(x, \mu) = \int_0^{L/|\mu|} e^{-\Sigma_t \rho} d\rho = \frac{1 - e^{-\Sigma_t L/|\mu|}}{\Sigma_t}$$

Now integrate over μ:

$$\int_0^1 J(x, \mu)\,d\mu = \frac{1}{\Sigma_t}\int_0^1 (1 - e^{-\Sigma_t L/\mu}) d\mu$$

The integrand `e^{-Σ_t L/μ}` is stiff near μ = 0 (exp of a large argument).

**Substitution**: `v = -ln(μ)`, so `μ = e^{-v}`, `dμ = -e^{-v} dv`:

$$\int_0^1 e^{-\Sigma_t L/\mu} d\mu = \int_\infty^0 e^{-\Sigma_t L e^v} (-e^{-v}) dv = \int_0^\infty e^{-\Sigma_t L e^v}\,e^{-v}\,dv$$

The integrand `e^{-v} · e^{-Σ_t L e^v}` is a product of standard Laguerre weight `e^{-v}` and a super-exponentially-decaying `e^{-Σ_t L e^v}`. Gauss-Laguerre on `v ∈ [0, ∞)` handles this naturally.

**Connection to E_1**: with one more substitution `u = Σ_t L e^v`, `du = Σ_t L e^v dv = u dv`, so `dv = du/u`:

$$\int_0^\infty e^{-\Sigma_t L e^v} e^{-v} dv = \int_{\Sigma_t L}^\infty e^{-u} \frac{e^{-v}}{u} du$$

With `e^{-v} = Σ_t L / u`:

$$= \int_{\Sigma_t L}^\infty \frac{e^{-u} \Sigma_t L}{u \cdot u} du = \Sigma_t L \int_{\Sigma_t L}^\infty \frac{e^{-u}}{u^2} du$$

Using integration by parts or the recurrence: `∫_a^∞ e^{-u}/u² du = (e^{-a}/a) - E_1(a)`:

$$= \Sigma_t L \cdot \left[\frac{e^{-\Sigma_t L}}{\Sigma_t L} - E_1(\Sigma_t L)\right] = e^{-\Sigma_t L} - \Sigma_t L \cdot E_1(\Sigma_t L)$$

Substituting back gives an expression combining exp and E_1 — which is PRECISELY the `E_2(Σ_t L)` that appears in slab problems.

### App C: Ki_1 → E_1/2 planar limit — sketch

The 2-D transverse Green's function `G_2D(ρ) = Ki_1(Σ_t ρ)/(2π ρ)` arose from z-integration of the 3-D point kernel. In the planar limit (`r_0 → ∞`, fixed L), rays in the cylinder stay close to the tangent plane. The z-direction becomes equivalent to a transverse direction in the slab.

Comparing the two kernels:
- Slab: `G_{\rm slab}(|x-x'|) = E_1(Σ_t |x-x'|)/2`
- Cylinder: `G_{\rm cyl}(ρ) = Ki_1(Σ_t ρ)/(2π ρ)`

At large `r_0`, rays from `r_i = r_0 + x` in direction `β` (with `x = const`) have transverse 2-D distance `|Δr| = 2r_0 |sin(β/2)| ≈ r_0 β` for small β. The cylinder formula `Ki_1(Σ_t · r_0 β)/(2π r_0 β)` is evaluated.

For small `β`, the integrand has `ρ = r_0 β`, and integrating `Ki_1/(2π r_0 β) · r_0 dβ = Ki_1(Σ_t r_0 β)/(2π β) dβ`. Substituting `μ = cos(angle from x-axis) ≈ 1 - β²/2`, and recognizing the definition `E_1(x) = 2 ∫_0^{π/2} sin(φ) Ki_1(x/sinφ) dφ / (... something)` ...

This derivation requires careful asymptotic matching and is left for the Phase-G.1 documentation commit. It's the kind of thing that SymPy could verify numerically: for a specific `r_0`, compute `Ki_1` integrated over the cylinder's polar angle and compare to the slab's `E_1/2` — the ratio should approach 1 as `r_0` → ∞.

### App D: Woodbury identity for rank-N BC update

For rank-N perturbation `K' = K + U V^T` where `U, V ∈ R^{N_r × N_modes}`:

$$(K + UV^T)^{-1} = K^{-1} - K^{-1} U (I + V^T K^{-1} U)^{-1} V^T K^{-1}$$

Cost: one `N_modes × N_modes` inverse (trivial) + two `N_modes × N_r` matrix-vector products. Total: O(N_modes² · N_r + N_modes³), vs rebuilding the full matrix at O(N_r³).

For the white-BC rank-N closure (Chapter 8), this lets us compute `(I - K_scatter)^{-1}` ONCE per eigenvalue iteration and update it with the BC rank-N perturbation efficiently.

---

## Part VI — Risk and cost-benefit analysis

### 6.1 Phase F (hollow-core)

| Aspect | Assessment |
|---|---|
| Effort | Medium (~1-2 days) |
| Risk of regression | Low — new code path (`inner_radius > 0`) orthogonal to solid case |
| Verification difficulty | Medium — need Sanchez/literature values for hollow-core geometries (Zotero MCP was flaky; literature-researcher may need to hunt) |
| User-visible benefit | HIGH — real fuel-pin geometries (central plenum, coolant channel) |
| Unlocks future work | Phase G planar-limit check |

**Verdict**: Strong candidate for first execution.

### 6.2 Phase G (slab unification)

| Aspect | Assessment |
|---|---|
| Effort | Medium (~1-2 days) |
| Risk of regression | MEDIUM — if polar-form slab is slower/less-precise than E_1 form, we risk complaining users |
| Verification difficulty | Easy — extensive Phase 4.1 slab tests to cross-check |
| User-visible benefit | Moderate — architectural cleanup, but users don't directly see it |
| Unlocks future work | Single-module simplification |

**Verdict**: Do AFTER Phase F so we can use the planar limit for a strong cross-validation. Gate on benchmark.

### 6.3 Phase H (τ-coordinate)

| Aspect | Assessment |
|---|---|
| Effort | Small-medium (~1 day) |
| Risk of regression | LOW — feature flag (`use_tau_coordinate`), default off initially |
| Verification difficulty | Easy — direct apples-to-apples comparison with ρ-integration on same problems |
| User-visible benefit | HIGH in specific cases — optically thick cells, multi-region, hollow cavities |
| Unlocks future work | Natural for Phase F hollow cavities (cavity becomes τ-no-op) |

**Verdict**: **Do FIRST**. Lowest risk, highest benefit-per-effort, helps F and G.

### 6.4 Recommended execution order

```
1. Phase H.1 (docs) + Phase H.2 (ρ(τ) walker) + Phase H.3 (tau-integration) + Phase H.4 (benchmark)
   Commit and verify.

2. Phase F.1–F.6 (hollow-core). Uses τ-coordinate internally where beneficial.
   Commit and verify.

3. Phase G.1 (coordinate-transforms docs extension) + Phase G.2 (slab in unified framework) + 
   Phase G.3 (regression) + Phase G.4 (planar-limit verification)
   Commit and verify.

4. Phase G.5 (conditional peierls_slab retirement, gate on benchmark)
   Only if G.3 benchmark shows the unified framework is competitive.
```

Total estimated effort: **3-5 days of focused work**, producing ~12-18 commits.

---

## Part VII — Next-session quick-start (if you're the future Claude picking this up)

1. Read this plan end-to-end. Pay special attention to Part V (Mathematical appendix) — these derivations are critical and would be costly to re-derive.

2. Verify the CP campaign (parallel plan) is complete: `git log --oneline | head -20` should show commits related to Phases A–E.

3. Verify `CurvilinearGeometry` in `orpheus/derivations/peierls_geometry.py` has the structure described in §0.1. If not, the parallel plan is not yet complete and this plan is premature.

4. Read `docs/theory/peierls_nystrom.rst` sections 0–10 (the foundation).

5. **File new issues** on GitHub before starting work:
   - Issue: "Add hollow-core support to `CurvilinearGeometry`" (→ Phase F)
   - Issue: "Slab in unified Peierls framework via polar form" (→ Phase G)
   - Issue: "τ-coordinate quadrature option for Peierls Nyström" (→ Phase H)
   - Issue: "Rich Sphinx documentation for coordinate transformations in transport" (→ content across phases)

6. **Execute in the order recommended in §6.4**: H first, then F, then G.

7. Each phase should **update `peierls_nystrom.rst`** with a new section. By the end, that document becomes the definitive reference for the unified polar-form Peierls Nyström architecture across all 1-D radial geometries, with rich coordinate-transform content.

8. After each commit, run the full Phase-4.2 regression suite + any Phase-A–E sphere/CP regression suites to ensure no behavior change for existing cases.

---

## Appendix A — Critical excerpts from the parallel plan (for redundancy)

(The parallel plan lives in `~/.claude/plans/` which is ephemeral on container rebuild. Key content duplicated here.)

### The unified Peierls equation

$$\Sigma_t(r)\,\varphi(r) \;=\; \frac{\Sigma_t(r)}{S_d}\!\int_{\Omega_d}\!d\Omega\!\int_0^{\rho_{\max}(r,\Omega)}\!\kappa_d(\Sigma_t\rho)\,q(r'(\rho,\Omega,r))\,d\rho + S_{\rm bc}(r)$$

### The 3-tier kernel hierarchy

```
Level 0: point            e^{-τ}/(4πR²)
Level 1: Peierls          E_1/2,  Ki_1/(2π),  e^{-τ}/(4π)     ← pointwise Nyström
Level 2: partial current  E_2,    Ki_2,        e^{-τ}          ← escape, surface source
Level 3: flat-source CP   E_3,    Ki_3,        e^{-τ}          ← region averaging
```

Differential identities: `E_n' = -E_{n-1}`, `Ki_n' = -Ki_{n-1}`.

### Current-rank-1 BC approximation error (before Phase I higher-rank)

| R [MFP] | k_eff error vs k_inf |
|---|---|
| 1 | 21% |
| 2 | 7% |
| 5 | 2% |
| 10 | 1% |

### Key code files at post-CP-campaign entry

```
orpheus/derivations/
    peierls_geometry.py       Unified CurvilinearGeometry
    peierls_cylinder.py       Thin facade
    peierls_sphere.py         Thin facade (post Phase A)
    peierls_slab.py           NOT unified (yet)
    cp_geometry.py            Unified flat-source CP (post Phase B)
    cp_slab.py                Thin facade (post Phase B)
    cp_cylinder.py            Thin facade (post Phase B)
    cp_sphere.py              Thin facade (post Phase B)
    _kernels.py               mpmath E_n, Ki_n; BickleyTables RETIRED post Phase B
```

---

**End of plan. ~900 lines. Preserves session-N's topological-unification and coordinate-transformation insights for session-N+2 or later.**

---

## Appendix B — Post-CP-campaign update (2026-04-18)

Added AFTER the CP campaign (Phases A, B, C) landed on `feature/mms-broad` (commits `435c0b3`…`a4098dc`). The original plan was written at end of session-N when the CP unification was still a plan; this appendix captures what the campaign actually did, where its predictions were right/wrong, and patterns that emerged. Future sessions picking up Phases F/G/H should read this before trusting the original plan's cost/risk numbers.

### B.1 Bickley retirement postmortem — prediction vs reality

The plan's §7.1 (risk assessment) predicted a cylinder `k_inf` drift "up to ~5e-5" from retiring BickleyTables. **Actual drift: up to ~4.3e-4** — one order of magnitude larger than predicted.

Root cause (discovered during Phase B.4): the prediction assumed only the derivation-side tabulation error (~1e-3 per Ki_3 evaluation) was in play. **It missed that the runtime solver had its own independent Ki_3 tabulation with a DIFFERENT bias** — `_build_ki_tables` in `orpheus/cp/solver.py` used `cumsum * dx`, an O(h) left-Riemann sum that gave ~8e-4 per-kernel error. The solver's bias and the derivation's bias had been partially cancelling; removing the derivation's ~1e-3 bias exposed the solver's ~8e-4 bias as naked solver-vs-derivation disagreement.

**Resolution**: Phase B.4 ALSO rewired the solver to consume `_ki3_mp` from `cp_geometry` (same Chebyshev interpolant as the derivation). Solver and derivation now share bit-identical kernel evaluations — an entire class of "kernel-split" bugs is structurally prevented.

**Lesson for Phase F/G/H**: whenever a plan retires or replaces a kernel/primitive on the derivation side, **audit the runtime solver for its own copy of that kernel**. The ORPHEUS architecture before Phase B.4 had the solver maintaining private copies of Ki_3 / chord quantities — each a potential site of cancelled bias. Phase C's `chord_half_lengths` lift + Phase B.4's `_ki3_mp` sharing are the anti-pattern remediations. Future phases should continue this: **solver depends on derivations, never the reverse, and never with duplicated primitives.**

### B.2 Sphere Peierls postmortem — Issue #100's motivating claim retracted

The plan's §0.3 cites "Sphere issue #100" as evidence that rank-1 white-BC closure fails STRUCTURALLY on the sphere (`k_eff ≈ 6.7` for a bare 1G 1-region case vs `k_inf = 1.5`). Phase A of the CP campaign reproduced the sphere Peierls reference through the unified `CurvilinearGeometry(kind="sphere-1d")` and found this claim was **wrong**.

The unified implementation gives (for bare sphere, `k_inf = 1.5`):

| R/MFP | k_eff (unified, white) | err vs k_inf |
|---|---|---|
| 1  | 1.096 | 27 %     |
| 2  | 1.391 | 7 %      |
| 5  | 1.490 | 0.7 %    |
| 10 | 1.496 | 0.3 %    |

— i.e. sphere rank-1 behaves essentially identically to cylinder rank-1. The old `k_eff ≈ 6.7` pathology was a **missing `R²` cell-surface divisor** (the pre-Phase-A code used `R` for both cylinder and sphere, undercounting the sphere by a factor of R). The `CurvilinearGeometry.rank1_surface_divisor(R)` method added in commit `9d03948` is the fix.

**Implication for Chapter 8 (higher-rank BC)**: the motivation from §8's "Issue #100 structural failure" needs updating. Rank-1 fails **gradually** at thin R in both cylinder and sphere — it's the same physics, not a sphere-specific catastrophe. The higher-rank closure is still the right architectural fix (unblocking `continuous_cases()` registration at thin R via Issue #105/N3), but it's an **accuracy upgrade**, not a **correctness fix**.

### B.3 New technique — scaled-kernel Chebyshev (missing from Part V)

The plan's §5.4 and §5.7 assume Chebyshev interpolation of Ki_3 will give ~1e-12 accuracy on [0, 50]. **It won't — plain Chebyshev on Ki_3 caps at ~3e-8 with 128 nodes** because the exp-decaying tail spans 22 orders of magnitude and a single polynomial can't resolve both the π/4 plateau at τ=0 and the exp(-50) tail.

**The fix used in Phase B.4**: interpolate the **scaled kernel** `f(τ) = exp(τ) · Ki_3(τ)` instead. The scaling converts the exp-decaying tail into a slowly-varying function that varies from 0.785 at τ=0 to 0.16 at τ=60. A degree-63 polynomial resolves `f` to ~2e-6 relative accuracy. Evaluate by multiplying by `exp(-τ)` after evaluation — preserves full accuracy on `Ki_3` itself (abs error ≤ 2e-6 · π/4 ≈ 2e-6 at τ=0, and ≤ 2e-6 · exp(-τ) elsewhere so relative error stays bounded).

**Code pattern** (see `cp_geometry.py::_ki3_scaled_cheb`):

```python
def func_scaled(tau):
    return np.array([float(ki_n_mp(3, t, 30)) * float(np.exp(t)) for t in tau])
poly = np.polynomial.Chebyshev.interpolate(func_scaled, deg=63, domain=[0, tau_max])
def _ki3(tau):
    return poly(tau) * np.exp(-tau)
```

This technique generalises to **any kernel with `exp(-τ)·g(τ)` structure**:
- Chapter 6's exp-stretched μ: same idea applied to the slab angular integrand.
- App A's τ-coordinate: the Gauss-Laguerre quadrature naturally absorbs `e^{-τ}`, which is the analytic equivalent of pre-scaling.
- Phase F hollow-core's ray through a cavity: the `e^{-τ_before_cavity}` factor is a constant across the cavity; factoring it out gives stable quadrature on the post-cavity segment.

**Add to Part V as App E**: "Scaled-kernel Chebyshev interpolation for exp-decaying functions." The derivation is trivial but the trick is not obvious — I (session-N+1) burned ~20 minutes trying plain Chebyshev at 128 nodes before pivoting. Preserve the lesson.

### B.4 Pattern — 2-commit refactor decomposition (relevant to F, G, H sequences)

Phase B.2 introduced `cp_geometry.py` and migrated three facades. The plan's §6 recommended two commits (one for new module, one for facade migration). **This turned out to be the right call** — `f1b869b` landed `cp_geometry.py` with 30 parity tests against unchanged legacy facades, proving bit-identity. `bf128d3` then migrated the three facades atomically. Bisection points are preserved; any Phase B.3/B.4 regression can be localised to one or the other commit.

**Generalisation for Phase F (hollow-core)**:
- Commit F.a: extend `CurvilinearGeometry` with `inner_radius` field, `rho_max` with inner-intersection, `optical_depth_along_ray` with cavity handling. **All existing tests pass** (default `inner_radius=0`).
- Commit F.b: add `compute_P_esc_inner` / `compute_G_bc_inner` + rank-2 BC closure; new tests exercise hollow geometry; solid geometry unchanged.
- Commit F.c: facade additions / `cp_geometry` hollow branch if applicable.

Each commit is individually revertible. The plan's F.1–F.7 sequence is already in this shape — keep it.

**For Phase H (τ-coordinate)**:
- Commit H.a: `optical_depth_along_ray_with_map` returning piecewise-linear ρ(τ) breakpoints. Testable in isolation; no callers change. (This is App A's Step 5 made executable.)
- Commit H.b: `use_tau_coordinate=True` flag in `build_volume_kernel`. Feature-flagged; `False` is default and unchanged.
- Commit H.c: benchmark test in `tests/derivations/test_peierls_tau_coordinate.py` — precision vs node count at R ∈ {1, 2, 5, 10, 20} MFP for both ρ and τ paths.
- Commit H.d (conditional, after benchmark): switch default based on `Σ_t * R > 5` heuristic (§5.5 prediction).

### B.5 Testing pitfall — nested mpmath.quad is quadratic-time

When writing verification tests for the kernel integration identities (Phase B.3 work on §12's inner/outer antiderivatives), my first draft used

```python
integral = mpmath.quad(
    lambda t: mpmath.quad(lambda u: ..., [0, mpmath.inf]),  # Ki_1 inner integral
    [a, b],
)
```

— i.e. mpmath adaptive quadrature of a function that itself required mpmath adaptive quadrature. **This ran for 45+ minutes per test before I killed it.** The outer quad's adaptive algorithm samples ~50 points, each of which triggers another ~50-point inner quad — O(2500) `exp(-t·cosh(u))/cosh(u)` evaluations per test case, with each evaluation already slow at 30 dps.

**The fix**: use `ki_n_mp(1, t, 20)` as the single-quad integrand (one quad, not two). `ki_n_mp` has its own internal quad but it's called once per outer sample:

```python
left = float(mpmath.quad(lambda x: ki_n_mp(1, float(x), 20), [a, b]))
right = float(ki_n_mp(2, a, 30)) - float(ki_n_mp(2, b, 30))
```

~0.5 s per test, same verification content.

**Add to Part V as App F**: "mpmath.quad composition pitfalls." Specifically, the rules:
- **Never** compose `mpmath.quad(mpmath.quad(...), ...)` inside a test.
- Numerical identities should be verified via **endpoint evaluation** of antiderivatives (right side) compared to **one-level** numerical integration (left side), not by performing the full chain of integrations on both sides.
- When in doubt, use closed-form equivalents (e.g. `mpmath.expint(n, x)` for E_n) instead of composed quads.

### B.6 MC cross-verification opportunity (under-developed in plan §5.3)

The plan's §5.3 notes that τ-coordinate is literally the Monte Carlo delta-tracking sampling variable, then lists three MC references. It does NOT develop the verification implication, which is genuinely valuable.

**Proposed new phase** (call it Phase H.5 or its own Phase J):

τ-coordinate Nyström and Woodcock-tracking MC are the **same integral** evaluated two ways:

- **Deterministic** (Nyström H): Gauss-Laguerre nodes at `{τ_k}_{k=1}^{N}`, weights `w_k`, deterministic evaluation.
- **Stochastic** (MC): sample `τ_i ~ Exp(1)` for each history, track until `τ = τ_max` or particle escapes, estimate via next-event estimator or track-length estimator.

Both produce an estimate of `∫_0^{τ_max} e^{-τ} · [integrand_in_τ_coords] dτ`. As MC sample count N_hist → ∞, MC estimate → Nyström estimate → exact, both governed by the same integral.

**This is a new verification band for ORPHEUS**: a cross-stochastic-deterministic L3-level comparison that doesn't currently exist. The `orpheus/derivations/mc.py` module already has the infrastructure; adding a Woodcock tracker that uses the same `CurvilinearGeometry` ray walker would make the two codes share their geometry primitives. Any bug in the shared primitives surfaces as an MC-Nyström disagreement that grows with N_hist (rather than shrinking as √N_hist as pure MC noise does).

**This is genuinely new V&V coverage**. Document as Part IV §24 "Monte Carlo connections" — expand the brief "§24.1 Delta-tracking ≡ τ-coordinate sampling" sub-section into a full discussion with worked example.

### B.7 Issue-state reconciliation (as of 2026-04-18)

The plan references "Issue #100" as structural failure motivation for Chapter 8. Actual issue state:

| Issue | State | Relationship to this plan |
|---|---|---|
| #94 (Bickley naming/retirement) | CLOSED by Phase B.4 | n/a |
| #100 (sphere rank-1 insufficient) | OPEN | motivating claim retracted per §B.2; issue should be re-scoped to "rank-1 accuracy upgrade" — same as #103 |
| #102 (chord_half_lengths lift) | CLOSED by Phase C | n/a |
| #103 (higher-rank white-BC = N1) | OPEN | Chapter 8's target |
| #104 (multi-group Peierls = N2) | OPEN | orthogonal to this plan |
| #105 (register continuous_cases = N3) | OPEN | blocked on #103 |
| #107 (flat-source CP unified = N6) | OPEN — closing in next merge | Phase B delivered |
| #108 (Phase 4.3 sphere) | CLOSED by Phase A | n/a |
| #55 (2D ray-tracing CP) | OPEN | **τ-coordinate is the natural primitive** — note in Phase H commit body |
| #56 (interface current) | OPEN | orthogonal |
| #57 (group loop vectorisation) | OPEN | orthogonal; Phase B.4 may have accidentally helped |

**Missing from the plan**: Issue #55 (2D ray-tracing CP) is directly enabled by Phase H — τ-coordinate is the natural geometry-independent primitive for any future 2D-ray code. The plan should reference #55 in Phase H's motivation and in Part IV §24 (this is the MC connection again — MC codes have used τ-space ray tracing for decades because 2D/3D geometry made ρ-space ray tracing impractical).

### B.8 Sphinx-build discipline (add to §0 prerequisites)

The original plan's §0 lists what to read before executing. **Add to that list**: "Run `sphinx-build -W docs docs/_build/html` before starting and verify clean. Run it after every commit and verify no new warnings."

The Phase B campaign dropped several stale-doc entries from the Nexus graph simply by rebuilding after the commits that modified the referenced symbols. The archivist also uses `-W` mode to detect duplicate-citation warnings (pre-existing, not to be fixed) vs NEW warnings (to be fixed). This discipline catches doc-code drift at the commit boundary rather than accumulating.

### B.9 Recommended execution update

The plan's §6.4 recommends: **H → F → G**. Based on what we now know:

- **H.1 (τ-coord docs)** is the cheapest single commit that preserves mathematical content (App A is ~300 lines of Sphinx). **Do this first regardless** — the math is valuable even without implementation.
- **Chapter 8 → Issue #103** (higher-rank BC) is higher physics value than H.2–H.4 (τ-coord implementation), because rank-1's 27% error at R=1 MFP is a concrete verification gap today; τ-coord is an efficiency win with less urgency. **Consider swapping to: H.1 docs → #103 implementation → H.2–H.4 implementation → F.**
- **G (slab unification)** remains low-priority. With Phase B.4 having established the "solver and derivation share primitives" pattern, slab is the last module that doesn't yet share this — but there's no user-visible benefit and the exp-stretched μ quadrature is unproven. Keep deferred unless a specific user need emerges.

**Net recommendation** for next session: H.1 (pure docs), then #103 implementation (physics), then decide H.2–H.4 vs F based on user direction.

---

**End of Appendix B.** Total additions: ~220 lines, written 2026-04-18 by the session that completed the CP campaign (commits `435c0b3` sphere scaffold through `a4098dc` BickleyTables doc sweep). Preserves patterns and postmortems that the original plan's author could not have known.
