# Cylinder Variant α Phase 1b (Multi-Region) — Verification Plan

**Status**: Pre-implementation. Tests-first (L7 per `.claude/lessons.md`).
The method-implementer dispatch downstream will transcribe these gates
into pytest while building `solve_greens_function_cylinder_mr`.

**Context**: closes ERR-026 / Issue #168 cylinder-MR gap by promoting
the cylinder Variant α Green's-function solver from Phase 1
(homogeneous, `solve_greens_function_cylinder_mg` at
[`greens_function_cylinder.py:439`](../../orpheus/derivations/continuous/trajectory_resolvent/greens_function_cylinder.py))
to Phase 1b (multi-region), mirroring the sphere MR pattern at
`solve_greens_function_sphere_mr`
([`greens_function.py:754`](../../orpheus/derivations/continuous/trajectory_resolvent/greens_function.py)).

---

## §1 — Scope and what's being verified

### New functionality (relative to Phase 1 homogeneous cylinder)

The Phase 1b extension adds:

1. **Piecewise chord arithmetic**: a single in-plane chord of length
   `L_2D,period = 2√(R²-b²)` (impact parameter `b = r·|sin φ_az|`) is
   subdivided into `K_chord ≤ K` annular segments by intersecting the
   chord with the K-1 interior interface circles. Per-segment
   arclengths Δs_k must satisfy Σ Δs_k = L_2D,period to machine
   precision.

2. **Piecewise optical-depth accumulation**: the homogeneous bounce-
   period attenuation `e^{-Σ_t L_period}` becomes `e^{-Σ_k Σ_t^{(k)}
   ΔL_period_k}` where `ΔL_period_k = Δs_k / s_in-plane`. Likewise for
   the first-leg attenuation `e^{-Σ_t L_0}`.

3. **Per-region scalar-flux output**: `phi_g_k(r)` resolved on the
   composite radial grid with `region_at_node[i]` tagging each radial
   node's region.

4. **Bounce-period operator on piecewise τ**: the homogeneous
   geometric series `1 / (1 - α e^{-Σ_t L_period})` generalises to
   `1 / (1 - α e^{-τ_period})` where `τ_period = Σ_k Σ_t^{(k)}
   ΔL_period_k`. Algebraic structure preserved; only the optical
   depth becomes piecewise.

5. **Per-region scattering, fission, χ**: full `(K, G, G)` SigS,
   `(K, G)` νΣ_f, `(K, G)` χ tensors. The source profile q(r) becomes
   discontinuous at interfaces.

### What is UNCHANGED (must still hold)

- Axial weighting `s_in-plane = √(1-μ_axial²)` and 3D path length
  rule `s_3D = s_2D / s_in-plane`.
- Full Gauss-Legendre quadrature on `μ_axial ∈ [-1,1]` and
  `φ_az ∈ [0, 2π)` (grazing endpoints never landed on exactly).
- α-bounce-sum closed form `ψ_surf = αB / (1 - α e^{-τ_period})`.
- Variant α grazing-ray cancellation (the 0/0 form at tangential
  `φ_az → ±π/2` for α=1).

### Mathematical claim hierarchy (per `vv-principles` §Hierarchical
### claim taxonomy)

| Gate | Layer claimed | Pillar |
| ---- | ------------- | ------ |
| 1 | Bit-identity (foundation invariant) | Reduction-to-homogeneous |
| 2 | Eigenvalue | Closed-form (WM-72 cross-pillar) |
| 3 | Eigenvalue | Closed-form (analytical k_∞) |
| 4 | Flux-shape (continuity invariant) | Foundation |
| 5 | Flux-shape / eigenvalue | Semi-analytical (literature, if exists) |
| 6 | Convergence-order | MMS-adjacent / refinement |
| 7 | Operator-algebra (optional) | Foundation |

---

## §2 — V&V level matrix

| # | Gate | Level | Pillar | Failure mode caught | Structural bug detected |
| - | ---- | ----- | ------ | ------------------- | ----------------------- |
| 1 | Homogeneous-limit reducibility | foundation | Reduction-to-homogeneous | Mode #3 (missing factor), #5 (wrong index) | Per-region loop drops a segment / overcounts; piecewise τ accumulator wrong sign |
| 2 | WM-72 vacuum-BC cross-check | L1 | Closed-form (singular-eigenfunction Fredholm) | Mode #1 (sign flip), #3 (missing factor) | Wrong impact-parameter formula; wrong axial-weight propagation through MR machinery |
| 3 | Single-region specular k_∞ | L1 | Closed-form (k_∞ = νΣ_f / Σ_a) | Mode #2 (variable swap), #6 (convention drift) | SigS transpose; per-region νΣ_f indexing |
| 4 | Interface φ-continuity | foundation | Reduction-invariant | Mode #5 (wrong index) | Off-by-one in `region_at_node`; spurious flux jump from miscounted segment |
| 5 | Literature MR benchmark | L1 (if available) | Semi-analytical | Mode #1, #3, #6 | Composite optical-depth path error invisible to Gates 1-4 |
| 6 | Quadrature convergence | L1 | Refinement | Mode #4 (wrong recurrence) | Per-segment quadrature recursion wrong; non-monotone convergence |
| 7 | Reciprocity (optional) | foundation | Operator-algebra invariant | Mode #2 (variable swap) | Operator asymmetry between forward/backward chord traversal |

### Decorators per gate

| # | `@pytest.mark.<level>` | `@pytest.mark.verifies(...)` | `@pytest.mark.catches(...)` |
| - | ---------------------- | ---------------------------- | --------------------------- |
| 1 | `foundation` | `peierls-greens-cylinder-mr-homogeneous-reduction` | — |
| 2 | `l1` | `peierls-greens-cylinder-mr-wm72-vacuum` | — |
| 3 | `l1` | `peierls-greens-cylinder-mr-kinf` | — |
| 4 | `foundation` | `peierls-greens-cylinder-mr-interface-continuity` | — |
| 5 | `l1` | `peierls-greens-cylinder-mr-literature-benchmark` (if exists) | — |
| 6 | `l1` | `peierls-greens-cylinder-mr-quadrature-convergence` | — |
| 7 | `foundation` | `peierls-greens-cylinder-mr-reciprocity` | — |

No `@catches("ERR-NNN")` markers — this plan is **forward
verification**, not bug-reproducer regression. The ERR-026 catalog
entry should reference these gates in its "How it would have been
caught" section, but no gate carries an explicit `catches` decorator
because no extant bug is being trip-wired.

---

## §3 — Gates in priority order

### Gate 1 (foundation) — Homogeneous-limit reducibility

**Purpose.** Cardinal Bit-Identity criterion 1 (`vv-principles`
§ "Bit-identity vs principled-equivalence"): an MR run with all K
regions assigned the same material **MUST reproduce
`solve_greens_function_cylinder_mg`** to within `(K × n_traj_quad
× ULP × condition_number)` — the FP-non-associativity floor for the
piecewise integration tree.

**Why this is the load-bearing gate.** Single-region homogeneous is
already verified (Sood `Ua-1-0-CY` at 8.5e-6 in
`test_peierls_greens_function_cylinder_xverif_sood2003.py`). If
Gate 1 passes, **all of that verification chain inherits** into the
MR code path. Every subsequent gate then tests only the
*piecewise-specific* machinery (interface arithmetic, segment
accumulation) — not the underlying Green's-function physics.

**Configuration.**

- **1G case**: σ_t = 1.0, σ_s = 0.5, νσ_f = 0.75 (fuel-A).
  - K = 3 with `radii = [0.4·R, 0.7·R, R]`, all same material →
    compare against `solve_greens_function_cylinder_mg` (R = 5.0
    cm) at `n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32`.
  - K = 5 with `radii = [0.2·R, 0.4·R, 0.6·R, 0.8·R, R]`, all same
    material → same comparison.
  - Rationale for K=3 AND K=5: exercises both an interior and a
    growing accumulation chain. Bug in segment ordering would
    show one passing and the other failing.
- **2G case**: get_xs("A", "2g") for the same K=3 configuration.
  Must reduce to `solve_greens_function_cylinder_mg` with same XS
  shape at machine precision.

**Tolerance.** `rtol=1e-10` on k_eff. Strict — bit-identity is
expected here because the MR code path with all-uniform material
should execute precisely the same arithmetic as the MG path
(per the V_α1 invariance under chord-length partitioning).

**The non-bit-identity escape clause.** If the K-region reduction
necessarily changes the FP reduction tree (e.g., the implementation
uses a per-segment `np.einsum` that orders summations differently),
the tolerance relaxes to `rtol = K * n_traj_quad * 1e-14` per
`vv-principles` §"Bit-identity" criterion 3 (FP-non-associativity
bound). MUST be principled — relaxation requires documenting:
(a) the named intermediate that changed, (b) the structurally-
independent reference at the converged value, (c) the drift bound.

**Anti-pattern flag.** Using a homogeneous 1G config makes Gate 1
blind to the SigS-transpose-on-MR-broadcast variant of Mode #6
(convention drift). The 2G sub-gate is mandatory.

### Gate 2 (L1) — WM-72 vacuum-BC cross-check

**Purpose.** Structurally-independent L1 reference. WM-72 (Westfall-
Metcalf 1972) uses **singular-eigenfunction Fredholm iteration with
Wiener-Hopf factorization** — no Bickley-Naylor, no bouncing-chord
trajectory integrals. The cylinder Variant α uses
**angle-resolved Green's-function integration along bouncing
characteristics**. These are mathematically distinct frameworks
(different structural ground per `vv-principles` § "Structural
independence"); agreement implies both are converging to the same
true value.

**Configuration.**

- 1G isotropic at `c = 1.30` (Sood `Ua-1-0-CY`).
- σ_t = 0.32640 cm⁻¹, σ_s = 0.248064, νσ_f = 0.176256.
- **Single-region MR shape** with `radii = [R_c]` where R_c is
  swept to find criticality (driven by α=0, vacuum BC):
  `radii=np.array([5.284935])`, then solve with α=0 → expect k_eff
  = 1.0.
- For the actual cross-check at α=0: solve `solve_singular_eigenfunction
  _cylinder_bare_critical(c=1.30, sigma_t=0.32640, n_grid=24)` to
  obtain `r_c_cm`. Then run `solve_greens_function_cylinder_mr` with
  `radii=[r_c_cm], alpha=0.0` and require `|k_eff - 1.0| ≤ 1e-5`.
- Quadrature: `n_r=24, n_mu_axial=16, n_phi_az=32, n_traj_quad=64`.

**Tolerance.** `|k_eff - 1.0| ≤ 1e-5` (matches the existing
homogeneous cross-check ceiling). The achievable precision is
≤ 1e-5 abs (research-grade); the WM-72 path itself reaches ~3e-7
which exceeds the 1e-5 brief target.

**Why ≥ Gate 1?** Even if Gate 1 passes bit-identity, the
single-region MR path might still mishandle the vacuum BC if the
piecewise-τ accumulator has a phase-of-segments dependency that
cancels for uniform material (Gate 1) but mis-routes through the
α=0 branch. Gate 2 forces α=0 (vacuum) at single-region MR — a
mode the homogeneous reduction does not exercise.

**Anti-pattern flag.** Do NOT compare against `solve_greens_function
_cylinder_mg(alpha=0.0)` — that would be procedurally independent
but structurally identical (same Variant α algorithm). The WM-72
cross is the structurally-independent pillar.

### Gate 3 (L1) — Single-region specular-BC k_∞ recovery

**Purpose.** Closed-form L1: α=1 single-region cylinder is
infinite-medium-equivalent because specular BC + isotropic infinite
cylinder = no leakage = k_∞.

**Configuration.**

- 2G case (the 1-group-degeneracy rule: 1G is degenerate for eigenvalue claims):
  use `analytical.derive_2g` Mixture A 2G XS or build a custom 2G
  matrix with asymmetric SigS (so SigS ≠ SigS^T — guards against
  Mode #6 convention drift). Suggested asymmetric 2G:

  ```
  Σ_t = [1.0, 1.5]
  SigS = [[0.4, 0.3],   # group 0 → group 1 downscatter 0.3
          [0.1, 0.8]]   # group 1 → group 0 upscatter 0.1
  νΣ_f = [0.6, 0.9]
  χ = [1.0, 0.0]  # all-fast emission
  ```

  Compute k_∞ analytically: `k_∞ = max eig(A⁻¹F)` where
  `A = diag(Σ_t) - SigS^T`, `F = χ ⊗ νΣ_f`. Use
  `orpheus.derivations.common.eigenvalue.kinf_homogeneous` or
  `analytical.derive_2g` for the reference value.
- Single-region MR: `radii=[5.0]`, α=1.0, all XS the 2G above.
- Quadrature: `n_r=16, n_mu_axial=12, n_phi_az=24, n_traj_quad=32`.

**Tolerance.** `|k_eff - k_∞| / k_∞ ≤ 5e-4` per brief. The
existing 1G/MG cylinder solver achieves machine precision at α=1
(closed cylinder is an exact k_∞ recovery — this is V_α1's
algebraic identity); the MR shim should preserve this, so a
tighter `rtol=1e-8` is the actual expectation. The 5e-4 is the
permissive ceiling; achieving ≤ 1e-6 is the target.

**Why an asymmetric SigS matters.** The L0-SN-009 (ERR-002)
catalog entry warns that symmetric scattering matrices make
`SigS = SigS^T` and the convention-drift bug invisible. Same
applies here: a 2G symmetric SigS in Gate 3 would hide an MR-only
broadcast that flips `SigS[k, :, :]` orientation.

**Why 2G not 1G.** Mandatory by the `vv-principles` 1-group-degeneracy rule
§ "1-group degeneracy". 1G's k_∞ = νΣ_f/Σ_a is shape-independent;
any errors in the spatial / angular discretization or in the
per-region tensor broadcast are invisible.

### Gate 4 (foundation) — Interface φ(r) continuity

**Purpose.** The scalar flux `φ(r) = ∫ψ(r,μ,φ_az)dμdφ_az` is
continuous across material interfaces — the angular flux ψ is
discontinuous (its derivative jumps because the source profile
jumps), but the integral over the angular domain averages out the
jump.

**Configuration.**

- 3-region 1G cylinder: `radii = [0.5·R, 1.0·R, R]`, R = 4.0 cm.
- Strong material contrast: σ_t = [2.0, 0.5, 1.5] (10× contrast
  across the inner/middle boundary).
- σ_s = [1.8, 0.3, 1.0], νσ_f = [0.1, 0.0, 0.05] (multiplying core,
  scattering moderator, lightly multiplying reflector).
- α = 1.0 (closed cylinder; vacuum sub-test optional with α=0).
- Quadrature: `n_r=36, n_mu_axial=16, n_phi_az=32, n_traj_quad=48`.
  Use a composite radial grid with refinement near interfaces if
  the prototype supports it.

**Test.** Extract `phi[i]` for the two radial nodes nearest each
interface (one on each side). Compute the interpolated values at
the exact interface radius `r_k` from the left side and from the
right side. Assert:

```
|φ(r_k⁻) - φ(r_k⁺)| ≤ 10 × quadrature_floor
```

where `quadrature_floor` is the residual error of the angular
quadrature against a smooth test field at the same n_mu_axial,
n_phi_az setting (empirically ~ 1e-5 for `n_mu_axial=16,
n_phi_az=32`). So the bound is ~ 1e-4 absolute on `|φ⁻ - φ⁺|`
relative to `max(φ⁻, φ⁺)`.

**Why this catches what Gate 1-3 miss.** Gates 1-3 are
spatially-aggregated (k_eff or k_∞) — they integrate over the full
domain and lose pointwise interface information. A bug in
`region_at_node` (off-by-one assigning the wrong region to a node
near an interface) is invisible to spatial integrals but produces a
visible jump in φ near the interface.

**Anti-pattern flag.** Do NOT use a symmetric material profile
(e.g., σ_t = [1, 2, 1] with reflective symmetry about r = R/2).
A sign flip in the segment-accumulator would produce a
**symmetric** error pattern that averages to zero in the
spatially-integrated norm. Pick an asymmetric profile (e.g.,
[2.0, 0.5, 1.5]) where each interface jump has a different
magnitude.

### Gate 5 (L1) — Published literature MR benchmark

**Status: TENTATIVE — REQUIRES LITERATURE SCAN.**

**Goal.** Find a published cylinder multi-region benchmark with
tabulated k_eff or flux values to provide a third structurally-
independent pillar (sphere has Garcia 2021; cylinder analogue is
the open question).

**Candidates from `scratch/literature/`:**

- **Sanchez-Ganapol 1983** — "Benchmark Values for monoenergetic
  neutron transport in one-dimensional cylindrical geometry with
  linearly anisotropic scattering." Bare-critical 1G P_1 anisotropic.
  **Very likely NOT multi-region** — abstract says
  "monoenergetic ... bare cylinder." Method-implementer MUST verify
  by reading §1-2 of the PDF. If MR is out of scope, Gate 5 reverts
  to "no published cylinder MR benchmark identified at this scope."
- **Garcia 2006** — "Collision Probabilities in r-θ-z Geometry."
  Probably escape probabilities not k_eff. Inspect for tabulated MR
  results.
- **Calvik 1965 / 1967** — "Collision Probabilities ... cylindrical."
  Same — CP, probably no MR k_eff tables.
- **Atkinson 1972** — "Numerical solution of Fredholm integral
  equations of the second kind with singular kernels." Methods paper,
  likely no MR benchmarks.
- **Mitsis 1963** — bare-critical configurations only.

**Probable conclusion.** No published cylinder MR benchmark of the
specific shape "concentric annular regions + isotropic scattering
+ specular/vacuum BC + tabulated k_eff or φ(r)" exists. The
verification chain TERMINATES at Gates 1-4 + 6, with the
structurally-independent anchor being Gate 2 (WM-72 single-region
reduction) rather than a direct MR cross-check.

**If a benchmark IS found** (method-implementer dispatch should
**explicitly search Sanchez 1977, Calvik 1965/1967, Garcia 2006,
and Sood-Foster-Parsons 2003 § cylinder section** — and dispatch
literature-researcher if uncertain):

- Configuration: transcribe XS, geometry, BC convention.
- Compute the convention map (Variant α uses scalar flux
  `φ = 2π ∫ψ dμ dφ_az`; Garcia-family uses `φ = ∫ψ dμ`; account
  for any 2π or 4π factors).
- Tolerance per the benchmark's claimed precision; typically 2-5%
  at non-interface points, 10-15% at interface-adjacent points
  (per the Garcia 2021 cubic-spline-smoothing pattern in the sphere
  MR cross-check).

**Deliverable for method-implementer.** Either ship Gate 5 as a
parametrised L1 test (if benchmark found) OR document in the test
file's module docstring:

> "No published cylinder MR benchmark of the appropriate scope was
> identified during the Phase 1b verification work. The chain of
> structural independence for cylinder MR Variant α terminates at:
>  Gate 2 (WM-72 single-region) + Gate 3 (analytical k_∞) +
>  Gate 4 (interface continuity invariant) + Gate 6 (quadrature
>  refinement). This is documented as a known V&V limitation; a
>  future search of Russian and Japanese journal archives, plus
>  the Sanchez–Pomraning textbook (1989), might surface a
>  shippable benchmark."

### Gate 6 (L1) — Quadrature refinement convergence

**Purpose.** All three quadrature axes (`n_traj_quad` for chord,
`n_phi_az` for azimuthal, `n_mu_axial` for axial) must exhibit
spectral-rate convergence on the smooth integrand away from grazing
rays. Failure of any one direction indicates a wrong-segment-
arithmetic bug that degrades the integrand smoothness.

**Configuration.**

- Use a "stress" configuration **different from any Gate 5
  candidate's exact configuration** (anti-pattern: testing
  against the same input twice — see §4). Suggested:
  - 3-region 2G: `radii = [0.3·R, 0.65·R, R]`, R = 3.0 cm.
  - Region 0 (core): get_xs("A", "2g") — multiplying.
  - Region 1 (mid): get_xs("D", "2g") — strong scatterer.
  - Region 2 (outer): get_xs("B", "2g") — moderator.
  - α = 0.7 (partial reflection — neither limit).

**Test.** For each of three quadrature axes, refine independently
holding the other two at high values. Required sequence per axis:

```
n_chord_quad ∈ {16, 32, 64, 128}    at n_phi_az=64, n_mu_axial=32
n_phi_az    ∈ {16, 32, 64, 128}    at n_chord_quad=64, n_mu_axial=32
n_mu_axial  ∈ {8, 16, 32, 64}      at n_chord_quad=64, n_phi_az=64
```

For each axis, assert:

1. Successive |Δk_eff| values are **monotone decreasing** (rates
   improving).
2. The ratio |Δk_eff_i+1| / |Δk_eff_i| ≤ 0.5 between adjacent doublings
   (super-algebraic — Gauss-Legendre on smooth integrand should give
   geometric convergence; 0.5 is the conservative bound).

**Tolerance.** The above are necessary conditions. Sufficient
condition is that the highest-resolution result agrees with a
cross-checked target (e.g., a Richardson-extrapolated value from
the same sequence). Failure of monotonicity in any axis signals a
non-smooth integrand contribution — typically a missing per-segment
piece or a wrong handling of segment boundaries.

**Anti-pattern flag.** Do NOT also test refinement convergence at
the Gate 5 benchmark's exact configuration — that would be
"verifying the code with the test's own input." Pick a different
config to keep the refinement test independent of the benchmark
agreement test.

**Cross-cutting hygiene per `numerical-bug-signatures` H4.**
Convergence rate is necessary, NEVER sufficient. This gate is
**not** an L1 correctness claim on its own — it is a **necessary
precondition** for the L1 claims in Gates 2, 3, 5 to be meaningful.
A solver that converges rapidly to the wrong asymptote is still
wrong. The L1 evidence comes from Gate 2 + Gate 3; Gate 6 only
proves the convergence is well-behaved.

### Gate 7 (foundation, OPTIONAL) — Reciprocity

**Purpose.** The Green's-function operator `G` should satisfy a
reciprocity relation: `⟨Gψ, φ⟩ = ⟨ψ, G^T φ⟩` for smooth ψ, φ.
For the cylinder operator this is the statement that swapping the
source-detector pair leaves the response unchanged.

**Status.** Defer to method-implementer. Check whether Phase 1
(`solve_greens_function_cylinder_mg`) already has a reciprocity
test; if YES extend to MR; if NO defer to a follow-on issue.

**Why optional.** Reciprocity in cylinder Variant α involves
swapping (r_source, r_detector) AND simultaneously reversing the
chord traversal direction. The Phase 1 V_α1 algebraic identity
already implies reciprocity at the closed-cylinder limit (via the
ψ_surf = q/Σ_t invariance under L_first); the MR extension does
not change the underlying symmetry. So Gate 7 is more of a
secondary invariant than a new claim.

---

## §4 — Anti-patterns to flag

### Per-gate anti-patterns

**Gate 1.**
- **Symmetric material profile** (uniform → no asymmetry to expose
  segment-ordering bugs). MITIGATED by the brief's "all K regions
  same material" being the **goal** of Gate 1 (it's the limit, not
  a test config). The accumulator test does need a non-trivial K
  sweep (K=3 AND K=5).
- **1G only**: hides Mode #6 SigS convention drift. MUST include 2G.

**Gate 2.**
- **Procedurally-independent comparison to MG path**: testing MR vs
  MG with α=0 catches the segment-arithmetic bug, BUT both use
  Variant α algorithm — only structurally independent if MG vs MR
  paths exercise different code. The WM-72 cross is the
  **structurally** independent pillar; the MG-vs-MR comparison is
  a useful Gate 1 sub-test but does not substitute for Gate 2.
- **Quadrature too coarse to reach 1e-5**: brief says n_traj_quad=64
  is needed; do not skimp.

**Gate 3.**
- **Symmetric SigS 2G** (e.g., SigS_01 = SigS_10): hides Mode #6.
  USE asymmetric SigS as listed above.
- **1G k_∞ test**: degenerate. The 2G case is mandatory.

**Gate 4.**
- **Symmetric material profile** (e.g., σ_t = [1, 2, 1]): sign-flip
  bug in segment-accumulator would average out spatially. USE
  asymmetric.
- **Tolerance too tight**: φ(r) at interface from a finite
  quadrature has an O(1/n_mu) angular-cancellation residual; the
  tolerance is `10 × angular_quadrature_floor`, NOT
  `numerical_precision`.

**Gate 5.**
- **Self-referencing benchmark**: do NOT use Variant α's own output
  as the "benchmark"; that's L4 cross-implementation, not L1.
- **Symmetric Cartesian-2D code as a "cylinder reference"**: 2D
  rectangular and 1D cylinder are different geometries; the
  conversion is not exact.

**Gate 6.**
- **Test refinement at Gate 5's exact configuration**: this is the
  "test the code with the test's input" tautology. Use a
  configurationally-different stress case.
- **Use of `np.allclose` at machine precision for convergence
  monotonicity**: monotonicity is a property of the **sequence**,
  not absolute precision; assert `|Δk_eff_i+1| < |Δk_eff_i|` and
  bound the rate, do not compare to machine ε.

### `vv-principles` Mode 7 (MMS simplification bias)

This plan does not propose MMS rows (it relies on closed-form
references). However, the design discipline of Mode 7 applies in
**how angular and radial discretisations are chosen for the
stress configuration** in Gate 6:

- **NEVER pick** a Gate 6 configuration with uniform radial
  density or with σ_t·R products that fall in a "comfortable
  middle" regime (1 ≤ σ_t·R ≤ 5). Pick one with significant
  contrast (σ_t·R range 0.5 to 30) so the chord-piercing
  arithmetic is genuinely stressed.
- **NEVER pick** a 2G case with symmetric down/up scattering.
  Inherited bias from textbook examples — use asymmetric SigS
  as in Gate 3 to expose per-region tensor-broadcast bugs.

### Architectural pattern risks (cylinder-specific structural
### concerns)

**1. Tangential grazing-ray pathology has K-region cancellation
structure.** In the homogeneous cylinder, V_α1 ensures
`ψ_surf = αB/(1 - α e^{-τ_period}) → q/Σ_t` is finite as the
tangential limit `φ_az → ±π/2` makes both numerator and
denominator vanish. In the MR case, B is a **piecewise** integral
and τ_period is a **piecewise sum**. The 0/0 cancellation must
happen across the piecewise structure — and the limit's
"finite value" is no longer q/Σ_t (single material), but a
**weighted average across regions** that depends on the chord's
region-traversal pattern. A naïve MR implementation that builds B
and τ_period separately and ratios them may produce numerical
0/0 with NaN where the homogeneous code produced clean
cancellation. **Gate 3 (single-region specular k_∞) is the
direct check** that the MR machinery preserves this cancellation
when K=1. If Gate 3 produces NaN or k_eff < 0, suspect the
grazing-ray handling first.

**2. Composite radial grid choice interacts with interface
sampling.** A single Gauss-Legendre on (0, R) for the radial grid
(as in the sphere MR) does NOT place nodes ON the interface
radii — nodes land in region interiors. This is fine for k_eff
(spatial integral over the full domain), but for Gate 4 (interface
continuity) it means the test must INTERPOLATE φ to the exact
interface radius from each side. The method-implementer must
either (a) expose `region_at_node` per-node and use that for
left-side / right-side splitting, or (b) ship a composite per-
region GL grid (each region has its own GL with endpoints near
the interface).

**3. Per-region χ broadcast in fission spectrum.** Sphere MR
allows `chi` to be `(n_regions, G)`, `(G,)`, or `None`. The
cylinder MR MUST replicate this exact API or risk a downstream
test parameter-passing bug. The cylinder MG has chi as `(G,)`
only — there's no precedent for `(n_regions, G)` in the cylinder
code path, so this is a new code path. Suggest method-implementer
copy the sphere MR's chi-broadcast logic verbatim.

**4. Per-region scattering matrix orientation.** The cylinder MG
uses `sigma_s[g_from, g_to]` (row-source). Sphere MR uses
`sigma_s[k, g_from, g_to]` (region-first). Cylinder MR MUST match
the sphere MR convention. The Gate 3 (asymmetric SigS) is the
specific test that catches a wrong convention here.

---

## §5 — File layout

### Primary test file

**`tests/derivations/test_peierls_greens_function_cylinder_mr.py`**

Parallel to `test_peierls_greens_function_mr.py` (sphere). Contains
Gates 1, 3, 4, 6 as foundation/l1 tests. Reuses the fixture
pattern from sphere MR (`@pytest.fixture(scope="module")` for
shared XS).

Suggested test function names (mirroring sphere MR):

- `test_mr_single_region_reduces_to_mg` (Gate 1, sub-test K=1
  bit-identity)
- `test_mr_K3_uniform_reduces_to_mg_1g` (Gate 1, K=3 1G)
- `test_mr_K5_uniform_reduces_to_mg_1g` (Gate 1, K=5 1G)
- `test_mr_K3_uniform_reduces_to_mg_2g` (Gate 1, 2G)
- `test_mr_single_region_kinf_2g_asymmetric_sigs` (Gate 3)
- `test_mr_interface_continuity_3region` (Gate 4)
- `test_mr_quadrature_convergence_chord` (Gate 6, n_traj_quad axis)
- `test_mr_quadrature_convergence_phi_az` (Gate 6, n_phi_az axis)
- `test_mr_quadrature_convergence_mu_axial` (Gate 6, n_mu_axial)

### Separate file for WM-72 cross-check

**`tests/derivations/test_peierls_greens_function_cylinder_mr_xverif.py`**

Gate 2 (WM-72 vacuum-BC cross-check) warrants its own file because:

- It has the longest runtime of any gate (WM-72 Brent root search +
  cylinder MR power iteration).
- Its purpose is the L1 structural-independence anchor — flagging
  it via filename clarifies that this is the L1 evidence, not a
  generic regression test.
- Convention: `_xverif_` suffix matches the established pattern
  (`test_peierls_greens_function_cylinder_xverif_sood2003.py`,
  `test_peierls_greens_function_xverif_ps1982.py`).

Single test function: `test_mr_single_region_vacuum_matches_wm72`.

### Optional separate file for Gate 5

**`tests/derivations/test_peierls_greens_function_cylinder_mr_literature.py`**
(only if a publishable benchmark is found per §3 Gate 5).

Parametrised test per benchmark configuration.

### Sphinx page

**`docs/theory/trajectory_resolvent.rst`** (existing) — extended
with a new section "Multi-region cylinder Variant α". See §6.

---

## §6 — Sphinx labels created

The following `:label:` strings need a corresponding `.. math::
:label:` block in `docs/theory/trajectory_resolvent.rst` (under
a new "Multi-region cylinder Variant α" section, parallel to the
existing sphere MR section anchored at `peierls-greens-mr-
trajectory-segments` line 1388 and `peierls-greens-mr-piecewise-
tau` line 1425).

**New labels (8):**

| Label | Math content | Test gate(s) |
| ----- | ------------ | ------------ |
| `peierls-greens-cylinder-mr-trajectory-segments` | Cylinder chord-piercing: chord at impact parameter b, intersections with K-1 interface circles at distances s_k = b·something → arclength partitioning Σ Δs_k = L_2D,period | Gates 1, 4, 6 |
| `peierls-greens-cylinder-mr-piecewise-tau` | Piecewise optical depth: `τ_period = Σ_k Σ_t^{(k)} · Δs_k / s_in-plane` | Gates 1, 2, 3, 4, 6 |
| `peierls-greens-cylinder-mr-bounce-sum-piecewise` | Bounce-sum operator on piecewise τ: `1 / (1 - α exp(-Σ_k τ_k))` | Gates 1, 2, 3 |
| `peierls-greens-cylinder-mr-homogeneous-reduction` | Reduction: when all Σ_t^{(k)} = Σ_t uniform, τ_period reduces to homogeneous Σ_t · L_period | Gate 1 |
| `peierls-greens-cylinder-mr-wm72-vacuum` | Statement: single-region MR vacuum BC critical radius matches WM-72 r_c to ≤ 1e-5 | Gate 2 |
| `peierls-greens-cylinder-mr-kinf` | Statement: single-region MR specular BC k_eff = k_∞ exactly (V_α1 cylinder, MR-extended) | Gate 3 |
| `peierls-greens-cylinder-mr-interface-continuity` | Statement: φ(r) continuous across interfaces despite ψ discontinuity | Gate 4 |
| `peierls-greens-cylinder-mr-quadrature-convergence` | Statement: each quadrature axis exhibits monotone convergence under doubling | Gate 6 |

**Optional (1, if Gate 7 ships):**

| `peierls-greens-cylinder-mr-reciprocity` | Operator reciprocity: ⟨Gψ, φ⟩ = ⟨ψ, G^T φ⟩ | Gate 7 |

**Optional (1, if Gate 5 ships):**

| `peierls-greens-cylinder-mr-literature-benchmark` | Statement: agreement with <reference> at <points> to <tolerance> | Gate 5 |

The method-implementer ships stub labels per `algebra-of-record`
§ "Sphinx stub vs rich narrative" — the archivist downstream expands
each label into the full derivation page.

---

## Summary for method-implementer brief

This plan ships **6 mandatory gates + 1 optional gate** for the
cylinder Variant α Phase 1b multi-region extension:

- **Gate 1 (foundation, bit-identity)**: K-region uniform-material
  → reduces to `solve_greens_function_cylinder_mg`. Multiple K
  values (3, 5) and groups (1G, 2G). The load-bearing inheritance
  gate.
- **Gate 2 (L1, structural-independence anchor)**: single-region
  MR vacuum BC matches WM-72 `solve_singular_eigenfunction_cylinder
  _bare_critical` for 1G isotropic c=1.30 to ≤ 1e-5. Separate file
  `..._mr_xverif.py`.
- **Gate 3 (L1, closed-form k_∞)**: single-region MR specular BC
  with asymmetric 2G SigS matches `kinf_homogeneous` to ≤ 5e-4
  (target ≤ 1e-6). 2G mandatory.
- **Gate 4 (foundation, interface invariant)**: 3-region asymmetric
  σ_t profile (10× contrast) → φ(r) continuous across interfaces
  to ~ 1e-4 abs.
- **Gate 5 (L1, tentative)**: published cylinder MR benchmark, if
  found via literature search. Likely terminates as "no benchmark
  identified."
- **Gate 6 (L1, refinement)**: each of three quadrature axes
  exhibits monotone convergence under independent refinement on a
  3-region 2G stress configuration.
- **Gate 7 (foundation, optional)**: reciprocity invariant.

The single highest-priority architectural risk is the **tangential
grazing-ray cancellation** in the MR-extended bounce-sum — Gate 3
is the direct probe. The single highest-leverage test is **Gate 1**
because it inherits the entire single-region verification chain
(8.5e-6 vs Sood) into the MR code path.
